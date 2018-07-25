library(optparse)
library(survey)
library(dplyr)
library(stringr)
library(broom)
# ---------------- Set Up ---------------- #
# set up options
optionList <- list(make_option(c("-d","--data"), action="store", type = "character", help="set data to be train or test"),
                   make_option(c("-p","--phenotype"), action="store", type = "character", help="set phenotype to be BMI or glucose"))

args <- parse_args(OptionParser(option_list=optionList))

# load data
# creates demographicVariables, ExposureDescription, NHData.test, NHData.train
load("assignment3.rdata")


# choose correct data set
if(args$data == "train"){
  data <- NHData.train
}else if(args$data == "test"){
  data <- NHData.test
}else{
  stop("Wrong Data Type")
}

if(args$phenotype == "BMI"){
  phen <- "BMXBMI"
}else if(args$phenotype == "glucose"){
  phen <- "LBXGLU"
}else{
  stop("Wrong phenotype")
}


data <- data[which(!is.na(data[,phen])),]

data <- data %>% mutate(sex = case_when( male == 1 ~ "male", 
                                         female == 1 ~ "female"),
                        race = case_when(white == 1 ~ "white",
                                         black == 1 ~ "black",
                                         mexican == 1 ~ "mexican",
                                         other_hispanic == 1 ~ "hispanic",
                                         other_eth == 1 ~ "others"))

dsn <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=T, data= data)
res <- data.frame(matrix(vector(), 0, 6))

for(exposure in ExposureList){
  if(phen == "BMXBMI"){
    f <- str_glue("scale(BMXBMI) ~ scale(log10({exposure} + 1e-07)) + RIDAGEYR + sex + race + INDFMPIR")
  }else{
    f <- str_glue("scale(log(LBXGLU)) ~ scale(log10({exposure} + 1e-07)) + RIDAGEYR + sex + race + INDFMPIR")
  }
  exp_name <- ExposureDescription %>% filter(var == exposure) %>% select(var_desc) %>% distinct(var_desc) %>% slice(1) %>% pull(var_desc)
  # Check there is match for two data
  if(length(table(data$sex[!is.na(data[,exposure])]))<=1 | length(table(data[,exposure])) <= 1){
    d <- data.frame(exposure, exp_name, phen, fit.estimate=NA,  fit.std.error=NA, fit.p.value=NA)
  }else{
    fit <- svyglm(f, data = data, dsn) %>% tidy() 
    target <- str_glue('scale(log10({exposure} + 1e-07))')
    if(target %in% fit$term){
      fit <- fit %>% filter(term == target)
      d <- data.frame(exposure, exp_name, phen, fit$estimate, fit$std.error, fit$p.value)
    }else{
      d <- data.frame(exposure, exp_name, phen, fit.estimate=NA,  fit.std.error=NA, fit.p.value=NA)
    }
  }
  res <- rbind(res, d)
}

colnames(res)<- c("Exposure ID", "Exposure Name", "Phenotype Name", "Estimate", "Standard Error", "P-value")

res$FDR <- p.adjust(res$`P-value`, method = "fdr")

output <- str_glue("{args$phenotype}_{args$data}.csv")

print(output)

write.csv(res,output, row.names = FALSE)