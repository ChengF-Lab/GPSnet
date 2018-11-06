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


### Code (for Matlab):
- Run the Raw_Module_Generation.m to generate the raw module, then run Module_Treat.m to treat the raw module. The two processes are very time-consuming, and we need to save the current results in file folder Data_mat/Raw_Module_Score and Data_mat/Raw_Module. The code need to be executed on computer cluster. 

- When accomplished these two processes, we can obtain about 60,000 raw module for each cancer type, and then run Module_Selected.m to obtain the final module for each cancer type. 

- Closet_Distance_Final.m is used to calculate the network proximity between the drug targets and the module gene.
