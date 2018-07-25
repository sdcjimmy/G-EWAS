
library(dplyr)
library(ggplot2)
library(qqman)
library(corrplot)

list.files("WTCCC/")

disease <- c('BD','CAD','CD','HT','RA','T1D','T2D')

write_all_data <- function(dis){
    d <- read.csv("hw1/results/T2D/22_filtered.csv")

    data <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(data) <- names(d)
    for(chr in 1:22){
        chr <- ifelse(chr < 10, paste0('0',chr),paste(chr))
        d <- read.csv(paste0("hw1/results/",dis,"/",chr,"_filtered.csv"))
        data <- rbind(data,d)
    }


    info <- read.delim(paste0("WTCCC/",dis,"/snps_info.tar.gz"), sep = '\t', header = FALSE, skip = 1)
    names(info) <- c("start","end","EGAV","WTCCCid","rsID")
    rsID <- data %>% select(WTCCCid) %>% left_join(info, by = "WTCCCid")
    data$WTCCCid <- rsID$rsID
    names(data)[1] <- "rsID"

    output_path <- paste0("hw1/results/all/",dis,".csv")
    write.csv(data,output_path)
    
    list(rsID, data)
}

#for(dis in disease){
#    write_all_data(dis)
#}

T2D <- write_all_data("T2D")
data <- T2D[[2]]
rsID <- T2D[[1]]



# Number of filtered SNPs for analysis
nrow(data)

head(data)

# Threshold
0.05/nrow(data)

significant_data <- data %>% filter(chi2_pvalue < 0.05/nrow(data))

# Manhanttan plot data preparation
man <- data %>% select(SNP = rsID, CHR = Chromosome, P = chi2_pvalue)
man$BP <- as.numeric(rsID$end)
man$CHR <- as.numeric(man$CHR)
man <- man[which(!is.na(man$BP)),]

manhattan(man)

qqplot(-log10(runif(nrow(data),0,1)), -log10(data$chi2_pvalue), xlim = c(0,6),
      xlab = "expectd -log10(pvalue) under null",
      ylab = "actual -log10(pvalue)", pch = 16, cex = .5 )

# Q-Q plot
qq(data$chi2_pvalue)

# Numer of SNPs exceeding bonferroni thresholds
p_adjust <- p.adjust(data$chi2_pvalue,method = "bonferroni")
sum(p_adjust <= 0.05)

# rs4506565
data %>% filter(rsID == 'rs4506565')

# Inflation factor
If <- qchisq(median(data$chi2_pvalue),1,lower.tail = F)/qchisq(0.5,1)

#install.packages("GenABEL")
library(GenABEL)
If <- estlambda(data$chi2_pvalue, method="median")$estimate

IFadusted.pvalue <- pchisq(qchisq(data$chi2_pvalue, 1, lower.tail = F)/If, df = 1, lower.tail = F)

# Manhanttan plot data preparation
man$P <- IFadusted.pvalue
manhattan(man)

p_adjust <- p.adjust(IFadusted.pvalue,method = "bonferroni")
sum(p_adjust <= 0.05)

#
all_data <- data.frame()
for(dis in disease){
    data <- read.csv(paste0("hw1/results/all/",dis,".csv"))
    data$disease <- dis
    if(nrow(all_data) == 0){
        all_data <- data
    }else{
        all_data <- rbind(all_data, data)
    }
}

qnlog10 <- function(p) {
        -log10(qunif(p[length(p):1]))
}
ggplot(all_data, aes(sample = -log10(chi2_pvalue))) + stat_qq(distribution = qnlog10) + geom_abline(col = 'red') + facet_wrap(~ disease)

all_data %>% group_by(disease) %>% filter(chi2_pvalue <= 1e-7) %>% summarise(SNPs_count = n()) 

all_data %>% group_by(disease) %>% filter(chi2_pvalue <= 1e-7) %>% summarise(SNPs_count = n()) %>% arrange(desc(SNPs_count))

all_data %>% filter(rsID == "rs6679677") %>%  filter(chi2_pvalue <= 1e-7)

(combn(disease,2))

cor_table <- matrix(rep(0,49),7,7)
colnames(cor_table) <- disease
row.names(cor_table) <- disease

all_data %>% select(rsID, chi2_pvalue)

rsID_list <- data.frame(rsID = c(), chi2_pvalue= c())
for(d in disease){
    r <- all_data %>% filter(disease == d) %>% dplyr::select(rsID, chi2_pvalue) %>% top_n(-200, chi2_pvalue)
    rsID_list <- rsID_list %>% bind_rows(r)
}
rsID_list <- rsID_list %>% distinct(rsID)
nrow(rsID_list)

for(i in 1:7){
    for(j in 1:i)
        if(i != j){
            d1 <- disease[i]
            d2 <- disease[j]
            d1 <- all_data %>% filter(disease == d1) %>% filter(rsID %in% rsID_list$rsID) %>% dplyr::select(chi2_pvalue) 
            d2 <- all_data %>% filter(disease == d2) %>% filter(rsID %in% rsID_list$rsID) %>% dplyr::select(chi2_pvalue)
            c <- cor(-log10(d1), -log10(d2), method = 'spearman')
            cor_table[i,j] <- c
            
        }
}

corrplot(cor_table, method = 'circle', type = 'lower')

rsID_list <- data.frame(rsID = c(), chi2_pvalue= c(), disease = c())
for(d in disease){
    r <- all_data %>% filter(disease == d) %>% dplyr::select(rsID, chi2_pvalue) %>% top_n(-200, chi2_pvalue)
    r$disease <- d
    rsID_list <- rsID_list %>% bind_rows(r)
}
nrow(rsID_list)

for(i in 1:7){
    for(j in 1:i)
        if(i != j){
            d1 <- disease[i]
            d2 <- disease[j]
            d1 <- rsID_list %>% filter(disease == d1) %>% dplyr::select(rsID)
            d2 <- rsID_list %>% filter(disease == d2) %>% dplyr::select(rsID)         
            cor_table[i,j] <- length(intersect(d1$rsID,d2$rsID))/200
            
        }
}

corrplot(cor_table, method = 'circle', type = 'lower')

head(all_data %>% filter(disease == disease[1]))

head(all_data %>% filter(disease == disease[5]))

all_data %>% filter(rsID == "rs6679677")

dim(all_data)

all_data %>% select(logp = chi2_pvalue) %>% mutate(logp = -log10(logp)) 

all_data %>% filter(chi2_pvalue != 0) %>% mutate(logp = -log10(chi2_pvalue)) %>% select(logp) %>% max

head(all_data)
