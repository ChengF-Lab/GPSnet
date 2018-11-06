function CD=Closet_Distance_Final(Cancer_Type);

load Data_mat/Gene_Distance_17W  %%%%% distance between genes
load Data_mat/Net_PPI
Net=Net(:,1:2);
Net(find(Net(:,1)-Net(:,2)==0),:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};
Genes=unique(Net(:));
Degree=Genes;
for i=1:length(Degree)
    Degree(i,2)=length(find(Net(:)==Degree(i)));
end
tic;
for II=1:272
    CD=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,Cancer_Type);
    toc
    disp(II)
end


function CD=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,Cancer_Type);
gene_1=Module_Selected(Cancer_Type,0.5,0.005);
load(['Data_mat/Gene_Drug_10uM']);
load Data_mat/Map_List
drug=Map_List{II,1};
[i1,i2]=ismember(drug,Drug_List);
ii=i2;
gene_2=Gene_Drug(find(Gene_Drug(:,2)==ii),[1 3]);  %%%% drug gene
gene_2=gene_2(:,1);
gene_1=intersect(gene_1,Genes);
gene_2=intersect(gene_2,Genes);
for i=1:1000
    g_1=drug_network_random_module_calculation(Degree,gene_1,i);
    g_2=drug_network_random_module_calculation(Degree,gene_2,i);
    CD(i,1)=CD_Cal(Genes,Distance,g_1,g_2);
end
CD(i+1,1)=CD_Cal(Genes,Distance,gene_1,gene_2);  %%%%% the first column
save(['Data_mat/Drug_Enrichment_',Cancer_Type,'/CD_',num2str(II)],'CD') 


function CD=CD_Cal(Gene,Distance,gene_1,gene_2);
[a,m1]=ismember(gene_1,Gene);
[a,m2]=ismember(gene_2,Gene);  %%% drug gene
d=Distance(m1,m2);
d=double(d);
CD=mean([min(d)]);


function Module_R=drug_network_random_module_calculation(Degree,Gene,II);
ctime=datestr(now,30);
tseed=str2num(ctime((end-5):end));
rand('seed',tseed*II);    %%%%% random seed reset
Degree=sortrows(Degree,-2);
Degree_Dis=unique(Degree(:,2));
for i=1:length(Degree_Dis)
    Degree_Dis(i,2)=length(find(Degree(:,2)==Degree_Dis(i,1)));
end
Degree_Dis=sortrows(Degree_Dis,-1);
DB=Degree_Dis(1,1)+1;  %%% degree bin
kk=Degree_Dis(1,2);
for i=2:length(Degree_Dis)
    kk=kk+Degree_Dis(i,2);
    if kk>200
        DB=[DB;Degree_Dis(i,1)];
        kk=0;
    end
end
for i=1:length(DB)-1
    x=find(Degree(:,2)<DB(i) & Degree(:,2)>=DB(i+1));
    NDB{i,1}=Degree(x,1);   %%%% node in each degree bin
end
Module_R=[];
[a,b]=ismember(Gene,Degree(:,1));
MBD=Degree(b(a),2);
MBD=[Gene,MBD];
for i=1:length(Gene)
    k=MBD(i,2);
    m=find(DB<=k,1,'first')-1;
    NB=NDB{m};
    Module_R=[Module_R;NB(ceil(rand*length(NB)))];
end