function [Enrichment,P_Value]=Module_Validation(Cancer_Type);
%%%%%%%%%%% check the enrichment of the module gene on several gene sets (SMGs, Associated gene and Survival gene)

alpha=0.5;
Cancer_Gene=Cancer_Module_Calculation(Cancer_Type,alpha,0.005);  %%%%% Cancer_Gene

load(['Data_mat/Driver_Gene/',Cancer_Type]);
load(['Data_mat/Associated_Gene/',Cancer_Type]);
load(['Data_mat/Survival_Gene/',Cancer_Type]);

load Data_mat/Net_PPI
Gene=unique(Net(:,1:2));

LL=100;

Enrichment=zeros(LL,3);
Enrichment(1,1)=length(intersect(Cancer_Gene,Driver_Gene));
Enrichment(1,2)=length(intersect(Cancer_Gene,Associated_Gene));
Enrichment(1,3)=length(intersect(Cancer_Gene,Survival_Gene));
for i=2:LL
    Random_Gene=Gene(randperm(length(Gene),length(Cancer_Gene)));
    Enrichment(i,1)=length(intersect(Random_Gene,Driver_Gene));
    Enrichment(i,2)=length(intersect(Random_Gene,Associated_Gene));
    Enrichment(i,3)=length(intersect(Random_Gene,Survival_Gene));
end

for i=1:3
    [mu,sigma] = normfit(Enrichment(2:end,i));
    P_Value(i)=1-normcdf(Enrichment(1,i),mu,sigma);
end

