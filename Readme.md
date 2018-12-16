## The data and Matlab code for paper "A Genome-wide Positioning Systems Network Algorithm for in silico Drug Repurposing".


### Data Description (Data_mat)
All data is saved as the .mat file because all the calculation is performed by Matlab. 

### Cancer_Specific_PPI/ folder:
We used the p-value<0.05 co-expressed PPI interactions to build the cancer type-specific PPI. 

### Mutation/ folder: 
- There are two columns for the data Mutation. The first column represents the Gene ID, and the second column represents the times of the number of the tumors with the corresponding mutated gene (mutation frequency). 

- Driver_Gene/ folder, Associated_Gene/ folder, Survival_Gene/ folder: The gene list of the driver gene, associated gene and survival gene for each cancer type.


### Net_PPI: 
Net_PPI is the human interactome (total ppi), and there are four columns. The first and second column are the gene id, the third column represents the co-expression correlation, and the fourth column represent the P-value of the correlation.

### Gene_Distance:
There are two parameters Distance and Genes. Distance is a 15137x15137 matrix, and the element represents the distance between the corresponding genes in Net_PPI, and Genes is the corresponding gene ID.

### Gene_Length:
Gene_Length has two column. The first column is the gene ID, and the second is the corresponding gene length. 

### Drug_Gene_10uM: 
It is the drug gene interaction. The first column of Gene_Drug is the gene ID, and the second column of Gene_Drug is the drug ID, and the corresponding drug name in the Drug_List.

### Map_List:
Map_List illustrates the drug and its ATC contents.

### Gene_Drug: 
It is the drug gene interaction. There are two parameters in this file (Gene_Drug and Drug_List). The first column of Gene_Drug is the gene ID, the second column of Gene_Drug is the Drug_ID, and the third colum of Gene_Drug is amplitude value. The corresponding drug name is in Drug_List, and there are 1309 drug in this dataset. 


### Code (for Matlab):
- "Raw_Module_Generation.m" to generate the raw module. 

- "Cancer_Module_Calculation.m" to obtain the final cancer module. 

- "Module_Validation.m" to check the enrichment of the cancer module gene on several 

- "Closet_Distance_ZScore.m" to calculate the network proximity between the drug targets and the cancer module gene. 

- "Drug_Gene_Set_Enrichment.m" to calculate the overlap enrichment between the drug targets and the cancer module gene.

### Tutorial: 
The code has been tested on the following systems: Linux Ubuntu 16.04 and Windows 7, with Matlab R2016b installed. 

Put all data in the "Data_mat/" folder first and create a "Raw_Module/" folder under "Data_mat/" folder. 

Run "Raw_Module_Generation.m" , and save raw module data in the "Data_mat/Raw_Module/". There are two parameters (Module and Score) in this intermediate result , where Module represents the gene set in each raw module and Score records the seed gene, final score and module size for the corresponding raw module. For we obtain about 60,000 raw module for each cancer type in our work, and this processes is very time-consuming. Than run ‘Module_Validation.m’ to check the cancer module from GPSnet. 

Run "Drug_Gene_Set_Enrichment.m" and "Closet_Distance_ZScore.m" to get the new indications for approved drugs using both gene set enrichment analysis and network proximity approaches respectively. 
