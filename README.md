# TG-GATES
Import and analysis of TG-GATES data for machine learning pipeline.

# Folders and Files
## Open-tggates_AllAttribute.zip
Compressed folder where single and repeat dose gene expression files (CEL) can be retrieved along with the atributes of each. Remember to download, unzip and save to local directory. 

## import_DEG_tggates.R
R code which details how to import TG-GATES gene expression files after initial download from database. Additionally, this file also details how to utilize the Affy package to remove batch effects from gene expressions and get differentially expressed genes through statistical approaches.

## NecDEG.csv
List of differentially expressed genes for necrosis causing chemical compounds. This list reduces the number of variables or features that will be analzyed for machine learning.