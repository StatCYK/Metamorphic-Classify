# Classifying Metamorphic versus Single-Fold Proteins with Statistical Learning and AlphaFold2


## Setup
### Basic requirements (Code was run and tested on LINUX machines.)
* python
* R packages "ranger", "pROC", "caret", "foreach", "doParallel"

### Install external packages

* SMICE: refer to [https://github.com/StatCYK/SMICE](https://github.com/StatCYK/SMICE) for installation

## Experiments

### Run SMICE
The config files containing the selected PDB IDs from the Atlas, CodNaS-Q, and the PDB databases are deposited in the folder `./configs`

See the [tutorial](https://github.com/StatCYK/SMICE) for the details of setting up the configuration file and conducting SMICE

### Extract Features
`./code/Feature_Extract.ipynb` provides the code for extracting features from the SMICE's predicted conformational ensembles and obtaining the training dataset as saved in the folder `./data/Metamorphic_classify.csv`

### Training Random Forest Classifier
run `./code/RF_train.R`
