
library(qdapRegex)
library(dplyr)
library(tibble)
library(glmnet)
library(ggplot2)

load("/n/groups/bmi704_spring18//assignment2//GSE40279_r2.Rdata")
load("/n/groups/bmi704_spring18//assignment2//GSE41169_r2.Rdata")

dim(gse40279.data)

get_age <- function(x){
    rm_between(x, 'age', 'y', extract=TRUE)[[1]]
}

train_age_title <- as.data.frame(gse40279.meta[,1])
train_age <- as.numeric(apply(train_age_title,1,function(x) rm_between(x, 'age','y', extract= TRUE)[[1]]))

test_age_title <- (gse41169.meta[,16])
test_age <- as.numeric(gsub("\\D","",test_age_title))

train_probe <- rownames(gse40279.data)
test_probe <- rownames(gse41169.data)
common_probe <- intersect(train_probe, test_probe)

train_data <- as.data.frame(gse40279.data) %>% rownames_to_column('probe') %>% filter(probe %in% common_probe) %>% column_to_rownames(var = 'probe')
test_data <- as.data.frame(gse41169.data) %>% rownames_to_column('probe') %>% filter(probe %in% common_probe) %>% column_to_rownames(var = 'probe')

dim(train_data)

summary(train_age)

sd(train_age)

summary(test_age)

sd(test_age)

dim(train_data)

dim(test_data)

impute_mean <- function(x){
    mean <- mean(x, na.rm = T)
    x[which(is.na(x))] <- mean
    x
}

train_data_imputed <- t(apply(train_data,1, impute_mean))
test_data_imputed <- t(apply(test_data,1, impute_mean))

set.seed(123)
alpha = 0.5
glmnet.Training.CV = cv.glmnet(t(train_data_imputed), train_age, nfolds=10,alpha=alpha,family="gaussian")
# The definition of the lambda parameter:
lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
# Fit the elastic net predictor to the training data
#glmnet.Training = glmnet(t(train_data_imputed), train_age, family="gaussian", alpha=0.5, nlambda=100) 

lambda.glmnet.Training

# Prediction on test data
pred <- predict(glmnet.Training.CV, s=lambda.glmnet.Training, newx=t(test_data_imputed))
# Prediction on training data
pred_train <- predict(glmnet.Training.CV, s=lambda.glmnet.Training, newx = t(train_data_imputed))

mean_train_age <- mean(train_age)
SS_total <- sum((train_age - mean_train_age)^2)
SS_res <- sum((train_age - pred_train)^2)
train_r_square <- 1 - (SS_res/SS_total)
train_r_square

cor(train_age, pred_train)

mean_test_age <- mean(test_age)
SS_total <- sum((test_age - mean_test_age)^2)
SS_res <- sum((test_age - pred)^2)
test_r_square <- 1 - (SS_res/SS_total)
test_r_square

cor(test_age, pred)

# Write the coefficient.csv data
coef <- coef(glmnet.Training.CV, s = lambda.glmnet.Training)
probe <- rownames(coef)
out_df <- data.frame(probe = probe, coef = as.numeric(coef))
write.table(out_df, "coefficeient.csv", col.names = F, row.names = F, sep = ',')

# Plot test correlation plot
qplot(pred, test_age, xlab = "Prediction Age of test samples", ylab = "True Age of testing samples")

# Plot training correlation plot
qplot(pred_train, train_age, xlab = "Prediction Age of training samples", ylab = "True Age of training samples")

test_accelarate <- pred - test_age
train_accelarate <- pred_train - train_age 

table(test_accelarate>0)

table(train_accelarate > 0)

test_sex <- (gse41169.meta[,'characteristics_ch1'])
test_sex <- ifelse(test_sex == 'gender: Male', 'M', 'F')

table(test_sex, test_accelarate>0)

chisq.test(table(test_sex, test_accelarate > 0))

train_sex <- gse40279.meta[,'characteristics_ch1.3']
train_sex <- ifelse(train_sex == "gender: M", 'M', 'F')

table(train_sex, train_accelarate > 0 )

chisq.test(table(train_sex, train_accelarate > 0 ))

library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()


rank_probes[-3]

selected_probes <- out_df[which(out_df[,2] != 0), ]
rank_probes_coef <- selected_probes[order(abs(selected_probes[,2]), decreasing = T),]
# Remove Intercept
rank_probes <- rank_probes_coef[-3,1]
genes <- getNearestGene(hm450[rank_probes])
genes <- unique(genes$nearestGeneSymbol)

genes

selected_train_data <- train_data_imputed[rank_probes,]
cor_matrix <- cor(t(selected_train_data))
pca <- prcomp(cor_matrix)

ggplot(as.data.frame(pca$x) ,aes(x = PC1, y = PC2)) + geom_point()

eigenvals <- pca$sdev^2
var_explained <- cumsum(eigenvals) / sum(eigenvals)
var_explained[2]
min(which(var_explained >= 0.9))

dim(pca$x)

sum(pca$sdev^2)

dim(selected_train_data)

cor_matrix <- cor(t(selected_train_data))

dim(cor_matrix)

dist <- dist((selected_train_data))
cluster <- hclust(dist)

head(rank_probes)

c <- cutree(cluster, 50)

head(c)

cluster1_probe <- names(which(c == 1))
cluster1_gene <- getNearestGene(hm450[cluster1_probe])
cluster1_gene
write.csv(cluster1_gene,"cluster1.csv")

cluster2_probe <- names(which(c == 2))
cluster2_gene <- getNearestGene(hm450[cluster2_probe])
cluster2_gene
write.csv(cluster2_gene,"cluster2.csv")

cluster3_probe <- names(which(c == 3))
cluster3_gene <- getNearestGene(hm450[cluster3_probe])
cluster3_gene
write.csv(cluster3_gene,"cluster3.csv")
