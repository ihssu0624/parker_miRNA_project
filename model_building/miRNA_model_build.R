###################################################################################
## Name: miRNA_model_building.R
## Goal: With specified miRNA and mRNA dataset, build an elastic net model and store model results
##       1)
## Author: Claire Su 
## Date Created: 02/15/2020
## Date Last Modified: 02/28/2020
## Notes:
###################################################################################
###################################################################################
## Loading Libraries
library(caret)
library(glmnet)
library(tidyverse)
###################################################################################
## Read in command arguments that specify which mRNA dataset to use and which miRNA to use as outcome
args = commandArgs(trailingOnly=TRUE)
x = args[1]
z = as.numeric(args[2])
print(x)
print(z)
###################################################################################
#  Read in normalized microRNA dataset 
#  Subset to the microRNA you want to build the model for specified by command line argurment
microrna_selected<-read.table("/nas/longleaf/home/isu/ms/miRNA_model_runs/seurat/data/microrna_selected.csv", fill = TRUE, header=TRUE,sep=',',stringsAsFactors = FALSE,quote="")
microrna_selected_new <- microrna_selected[,-1]
rownames(microrna_selected_new) <- microrna_selected[,1]
microrna_list<-rownames(microrna_selected_new)

# Create a directory to store output
model_results_path = paste0("/nas/longleaf/home/isu/ms/miRNA_model_runs/seurat/most_variable/output/", x,"/")
dir.create(model_results_path, showWarnings = TRUE, recursive = TRUE)
print(model_results_path)

# Read in mRNA dataset to build model on, specified by command line argument
mrna_selected_new<-readRDS(file=paste0("/nas/longleaf/home/isu/ms/miRNA_model_runs/seurat/data/mrna_selected_new_",x,".RDS"))
print(paste0("successfully loaded mrna filtered data, dimension: ", dim(mrna_selected_new)))

################################################
##### Model Building
################################################
# Set up 5 fold cross validation
control <- trainControl(method="repeatedcv", number=5, repeats=1)
seed <- 100
set.seed(seed)
print("set up training param")
print("the microrna is:")
print(microrna_list[z])

# Set up miRNA as outcome and transposed mRNA dataset as predictors
# Cell as rows and genes (miRNA/mRNA) as column
y=t((microrna_selected_new)[z,])
microrna_name<-microrna_list[z]
t_mrna<-t(mrna_selected_new)

# Split into training and testing
# Create data partitions for the data: input outcome vector, 1 partion, 70% in training 30% in test
inTraining <- createDataPartition(as.vector(y),  times = 1,p = .7, list = FALSE)
print("created partitions")
training_x <- t_mrna[ inTraining,]
testing_x  <- t_mrna[-inTraining,]
training_y <- y[ inTraining,]
testing_y  <- y[-inTraining,]
print("ds ready for modelin building")

# Center and Scale predictors before building elastic net model
preProcValues <- preProcess(training_x, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, training_x)

# Train elastic net model
model_glmnet <- train(x =trainTransformed, 
                      y =as.numeric(training_y), 
                      metric='Rsquared',
                      method = "glmnet",
                      trControl = control)

# Store model results but do not store datasets
new_model_glmnet <- model_glmnet
new_model_glmnet$trainingData<-NULL
new_model_glmnet$finalModel$call$x<-NULL

# Make prediction using the model built onto the testing dataset
preProcValues_test <- preProcess(testing_x, method = c("center", "scale"))
testTransformed <- predict(preProcValues_test, testing_x)

pred_micro <- predict(new_model_glmnet ,testTransformed )

# Set up dataframe to store correlation result from testing 
correlation_result <- data.frame(miRNA=microrna_list[z], correlation = NA)
correlation_result$correlation=cor(pred_micro, testing_y, method = "pearson")

## Save model summary
saveRDS(new_model_glmnet,paste(model_results_path,microrna_name, ".rds",sep = ""))

## Save correlation results
saveRDS(correlation_result, paste0(model_results_path,"correlation_results_",z,".rds"))

print(x)
print(microrna_name)
print("done!")

