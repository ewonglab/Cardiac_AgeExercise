# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Train and apply the cardiac transcriptomic age clock.
# Steps: 1) Load reference and study expression profiles; 2) fit the regression model; 3) predict age and export validation figures.
# =======================================
# Source reformatted from snRNA_AgeEx/Clock_prediction.ipynb.
# Set the variables in config/.env before running this script.
source(file.path("config", "paths.R"))

library(data.table)

meta= fread(publication_path("HEART_SNRNA_ROOT", "data/2024-04-03.Clock/GSE132040_MACA_Bulk_metadata.csv.gz"))
meta$colnames = paste0(meta$`Sample name` , '.gencode.vM19')  #match with the count colname
head(meta)

cnt = fread(publication_path("HEART_SNRNA_ROOT", "data/2024-04-03.Clock/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz"))

head(cnt)

keep = c('gene',meta$colnames[which(sapply(strsplit(meta$`source name`,"_"), `[`, 1) == 'Heart')])
keep

cnt = as.data.frame(cnt)
count = cnt[,keep]
rownames(count) = count$gene
count = count[,-1]

sample = meta[which(sapply(strsplit(meta$`source name`,"_"), `[`, 1) == 'Heart'),]
head(sample)

sample$age_month = as.numeric(sample$`characteristics: age`)
sample$sex = sample$`characteristics: sex`

sample$plate = sapply(strsplit(sample$`Sample name`,"_"), `[`, 3)
sample$letter = substring(sample$`Sample name`, 1, 1)

sample$condition = ifelse(sample$age_month <12, 'YS', 'MAS')
  
table(sample$plate)

sample$experiment = 'Tabula Muris Senis'

sample = sample[,c('source name', 'condition' , 'age_month', 'sex', 'experiment', 'plate' )]
colnames(sample)[1] = 'ID'
sample$sex = toupper(sample$sex) #to upper case
head(sample)

rm(cnt)

ex = read.csv(publication_path("HEART_SNRNA_ROOT", "data/2024-04-03.Clock/All_RNA_Pseudobulk_by_sample_Raw.csv"))
rownames(ex) = ex$X
ex=ex[,-1]
head(ex)

ex_meta = data.frame(ID= colnames(ex))
ex_meta$condition = c(rep(c('MAE', 'MAS', 'YE', 'YS'), c(2,4, 4,3)))
ex_meta$age_month = ifelse(ex_meta$condition %in% c('MAE', 'MAS'), 19, 10)
ex_meta$sex = sapply(strsplit(ex_meta$ID,"_"), `[`, 3)
ex_meta$experiment = 'ExerciseHeart'
ex_meta$plate = paste0('Run',sapply(strsplit(ex_meta$ID,"_"), `[`, 2))
ex_meta

#combine with the tabula muris senis meta file
allsample = rbind(sample, ex_meta)
head(allsample)

keep = intersect(rownames(count), rownames(ex))
allcount = cbind(count[keep,], ex[keep,])

head(allcount)

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))

  
  # Model performance metrics
data.frame(
  RMSE = RMSE,
  Rsquare = R_square
)
  
}

library(limma)

design <- model.matrix(~age_month + sex+experiment, data=allsample)
design

library(edgeR)

x.f <- allcount[which(apply(allcount, 1, function(xx, mc){length(which(xx>=mc))>=1}, mc=10)),]
message("Filter out ", dim(allcount)[1]-dim(x.f)[1], " genes")

y <- DGEList(counts=x.f, group=allsample$age_month, genes=rownames(x.f))
y <- calcNormFactors(y, method="TMM")

x.norm <- cpm(y, normalized.lib.sizes=TRUE)

treatment.design <- design[,c(1,2)]
dim(treatment.design)

batch.design <- as.data.frame(design[,-c(1,2)])
dim(batch.design)

corrected_vst_df <- removeBatchEffect(x.norm,design=treatment.design,covariates=batch.design)
corrected_vst_df[1:5,1:5]

library(ggfortify)

pca_res <- prcomp(corrected_vst_df, scale. = TRUE)

autoplot(pca_res, data=allsample, colour='age_month')

autoplot(pca_res, data=allsample, colour='sex')

autoplot(pca_res, data=allsample, colour='plate')

autoplot(pca_res, data=allsample, colour='experiment')

dat=t(corrected_vst_df)

set.seed(100)

index = sample(1:nrow(sample), 0.7*nrow(sample))  #only use the tabula muris to train the model
index

train = dat[index,] # Create the training data  row as mice and column as gene name
y_train=allsample[index,]

test = dat[-index,] # Create the test data row as mice and column as gene name
y_test=allsample[-index,]

dim(train)

library(glmnet)

lambdas <- 10^seq(2, -3, by = -.1)
set.seed(10)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(train, y_train$age_month, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# Best 
lambda_best <- lasso_reg$lambda.min 
lambda_best 
