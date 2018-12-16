function P_Value=Drug_Gene_Set_Enrichment(Cancer_Type);


alpha=0.5;
Cancer_Gene=Cancer_Module_Calculation(Cancer_Type,alpha,0.005);  %%%%% Cancer_Gene
load Data_mat/Gene_Drug
load Data_mat/Map_List
load Data_mat/Net_PPI
Gene=unique(Net(:,1:2));

Z_Score=zeros(length(Map_List),1);
P_Value=zeros(length(Map_List),1);
for i=1:length(Map_List)
    drug=Map_List(i,1);
    [i1,i2]=ismember(drug,Drug_List);
    x=Gene_Drug(find(Gene_Drug(:,2)==i2),[1 3]);
    x=x(find(abs(x(:,2))>0.67),:);
    LL=10000;
    Enrichment(1,1)=length(intersect(Cancer_Gene,x(:,1)));
    for j=2:LL
        Random_Gene=Gene(randperm(length(Gene),length(Cancer_Gene)));
        Enrichment(j,1)=length(intersect(Random_Gene,x(:,1)));
    end
    Z_Score(i,1)=(Enrichment(1,1)-mean(Enrichment(2:end)))/std(Enrichment(2:end));
    [mu,sigma]=normfit(Enrichment(2:end));
    P_Value(i,1)=1-normcdf(Enrichment(1,1),mu,sigma);    
end
