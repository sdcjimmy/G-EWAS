# Load Library
library(dplyr)
library(ggplot2)
library(optparse)

# Set options for bash
option_list <- list(
  make_option(c("-c","--chromosome"),action = "store",type = "character", help = "Number of chromosome"),
  make_option(c("-d","--disease"),action = "store", type = "character", help = "Disease type")
)

opt = parse_args(OptionParser(option_list=option_list))


# Load files
chr <- opt$c
dis <- opt$d

disease_path <- paste0("~/bmi704/WTCCC/",dis, "/Affx_gt_",dis,"_Chiamo_", chr,".tped.gz")
control1_path <- paste0("~/bmi704/WTCCC/58C/Affx_gt_58C_Chiamo_",chr,".tped.gz")
control2_path <- paste0("~/bmi704/WTCCC/NBS/Affx_gt_NBS_Chiamo_",chr,".tped.gz")

print("loading data ...")
print(control1_path)
control1 <- read.delim(control1_path, sep = '\t', header = FALSE)
control2 <- read.delim(control2_path, sep = '\t', header = FALSE)

control <- cbind(control1, control2[,-c(1:4)])
disease <- read.delim(disease_path, sep = '\t', header = FALSE)

output_path <- paste0("~/bmi704/hw1/results/", dis)


# Define function for SNPs processing
get_disease_allel_summary <- function(x){
    t <- table(x)
    freq <- c()
    
    for(i in 1:length(t)){
        allele <- unlist(strsplit(names(t)[i]," "))
        freq[allele[1]] <- ifelse(allele[1] %in% names(freq),freq[allele[1]] + t[i], t[i])
        freq[allele[2]] <- ifelse(allele[2] %in% names(freq),freq[allele[2]] + t[i], t[i])
     }
    
    freq <- sort(freq,decreasing = TRUE)
    
    if(length(freq) == 1){
        return (c(freq , 0, names(freq), NA))
    }else{
        return (c(freq[1],freq[2], names(freq)[1],names(freq)[2]))
   }  
}

get_control_allele_summary <- function(x){
    t <- table(x)
    freq <- c()
    homo <- c()
    heter <- 0
   
    
    for(i in 1:length(t)){
        allele <- unlist(strsplit(names(t)[i]," "))
        freq[allele[1]] <- ifelse(allele[1] %in% names(freq),freq[allele[1]] + t[i], t[i])
        freq[allele[2]] <- ifelse(allele[2] %in% names(freq),freq[allele[2]] + t[i], t[i])
        
        if(allele[1] != allele[2]){
            heter <- t[i]
        }else{
            homo <- append(homo, t[i])
        }
    }
    
    freq <- sort(freq,decreasing = TRUE)
    homo_maj <- max(homo)
    homo_min <- ifelse(length(homo)==1,0,min(homo))
    
    if(length(freq) == 1){
        return (c(freq , 0, names(freq), NA,homo_maj, 0 , 0))
    }else{
        return (c(freq[1],freq[2], names(freq)[1],names(freq)[2],homo_maj, heter,homo_min))
   }  
}

get_chi2_pvalue <- function(major_disease, minor_disease, major_control, minor_control){
    if(minor_disease <= 5 | minor_control <= 5 | major_disease <= 5 | major_control <= 5){
        p <- fisher.test(matrix(c(major_disease, minor_disease,major_control, minor_control),2,2))$p.value
    }else{
        p <- chisq.test(matrix(c(major_disease, minor_disease,major_control, minor_control),2,2))$p.value
    }
    p
}

get_HWD_pvalue <- function(Heter_major_count, Homo_count, Heter_minor_count, MAF_control){
    q <- MAF_control
    p <- 1-q
    obs <- c(Heter_major_count, Homo_count, Heter_minor_count)
    exp <- c(p*p, 2*p*q, q*q) * sum(obs)
    chi <- sum((obs - exp)^2/exp)
    p <- pchisq(chi, df = 1, lower.tail = FALSE)
    p
}


### Start processing the data

# get allele and HWE information for control group
control_freq <- as_data_frame(t(apply(control[,-c(1:4)], 1, get_control_allele_summary))) %>% 
    rename(Major_Allele_Control = V3, 
           Minor_Allele_Control = V4, 
           Major_Allele_count_Control = V1, 
           Minor_Allele_count_Control = V2, 
           Heter_major_count = V5, 
           Homo_count = V6, 
           Heter_minor_count = V7) %>%
    mutate(Major_Allele_count_Control = as.numeric(Major_Allele_count_Control), 
           Minor_Allele_count_Control = as.numeric(Minor_Allele_count_Control),
           Heter_major_count = as.numeric(Heter_major_count),
           Homo_count = as.numeric(Homo_count),
           Heter_minor_count = as.numeric(Heter_minor_count))

# get allele information for disease group
disease_freq <- as_data_frame(t(apply(disease[,-c(1:4)], 1, get_disease_allel_summary))) %>% 
    rename(Major_Allele_Disease = V3, 
           Minor_Allele_Disease = V4, 
           Major_Allele_count_Disease = V1, 
           Minor_Allele_count_Disease = V2) %>%
    mutate(Major_Allele_count_Disease = as.numeric(Major_Allele_count_Disease), Minor_Allele_count_Disease = as.numeric(Minor_Allele_count_Disease))

# The number of group numbers
n_control <- ncol(control) - 4
n_disease <- ncol(disease) - 4

# Output file format and process
out <- control[,c(2,1)]
names(out) <- c("WTCCCid","Chromosome")

output <- out %>% bind_cols(control_freq) %>% bind_cols(disease_freq) %>%
    mutate(MAF_disease = ifelse(Minor_Allele_Control == Major_Allele_Disease, Major_Allele_count_Disease/(n_disease*2), Minor_Allele_count_Disease/(n_disease*2))) %>% 
    mutate(MAF_control = Minor_Allele_count_Control/(n_control*2)) %>%
    mutate(odds_ratio = (Major_Allele_count_Disease * Minor_Allele_count_Control) / (Major_Allele_count_Control * Minor_Allele_count_Disease)) %>%
    rowwise() %>%
    mutate(chi2_pvalue = get_chi2_pvalue(Major_Allele_count_Disease, Minor_Allele_count_Disease, Major_Allele_count_Control, Minor_Allele_count_Control)) %>%
    mutate(HWD_pvalue = get_HWD_pvalue(Heter_major_count, Homo_count, Heter_minor_count, MAF_control)) %>%
    select(WTCCCid, Chromosome, Minor_Allele_Control, Major_Allele_Control, MAF_disease, MAF_control, odds_ratio, chi2_pvalue, HWD_pvalue) 



# Filtered  1. MAF in control less than 0.01 2. Deviate from HWE 3. Incorrect chromosome type
output_filtered <- output %>% 
    filter(MAF_control >= 0.01) %>%
    filter(HWD_pvalue >= 0.05) %>%
    filter(Chromosome == sub("^0+", "", chr))


# Output file

if(!dir.exists(output_path)){
    dir.create(output_path)
}

write.csv(output, file = paste0(output_path,"/",chr,"_all.csv"),row.names = FALSE)
write.csv(output_filtered, file = paste0(output_path,"/",chr,"_filtered.csv"), row.names = FALSE)