# parker_miRNA_project

This is a working progress of Claire Su's master's paper project for spring, 2020.

The goal of this project is to build elastic net models trained on bulk mRNA data to predict miRNA data in both bulk and single cell setting.

Organization of files:

1. bulk_RNA_processing: Contains python files that normalize bulk mRNA and miRNA data from TCGA, and matches them based on cell bar code.

2. sc_RNA_processing: Contains R script that processes single cell mRNA data to select the most variable genes used for predictors in the elastic net model.

3. model_building: Contains bash scripts that feeds in numbers that specify which miRNA to build the model for and which mRNA dataset to use. It also contains the R script that builds the elastic net model and performs 5-fold cross validation.