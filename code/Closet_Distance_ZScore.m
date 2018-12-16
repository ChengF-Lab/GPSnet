function Z_Score=Closet_Distance_ZScore(Cancer_Type,alpha);

load Data_mat/Gene_Distance  %%%%% distance between genes
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
if nargin < 2
    alpha=0.5;
end
for II=1:272
    Z_Score(II,1)=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,Cancer_Type,alpha);
    toc
    disp(II)
end

function Z_Score=S_AB_Cal_Block(II,Net,Degree,Genes,Distance,Cancer_Type,alpha);
gene_1=Cancer_Module_Calculation(Cancer_Type,alpha,0.005);  %%%%% Cancer_Gene
load(['Data_mat/Gene_Drug_10uM']);
load Data_mat/Map_List
drug=Map_List{II,1};
[i1,i2]=ismember(drug,Drug_List);
ii=i2;
gene_2=Gene_Drug(find(Gene_Drug(:,2)==ii),1);  %%%% drug gene
gene_1=intersect(gene_1,Genes);
gene_2=intersect(gene_2,Genes);
if isempty(gene_1) | isempty(gene_2)
    Z_Score=nan;
else
    for i=1:1000
        g_1=drug_network_random_module_calculation(Degree,gene_1,i);
        g_2=drug_network_random_module_calculation(Degree,gene_2,i);
        CD(i,1)=CD_Cal(Genes,Distance,g_1,g_2);
    end
    CD(i+1,1)=CD_Cal(Genes,Distance,gene_1,gene_2);  %%%%% the first column
    Z_Score=(CD(end)-mean(CD(1:end-1)))/std(CD(1:end-1));
end


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


%%%%%%% calculate the largest component of the PPI network
function [LG,L]=largest_component(G);
k=0;
if isempty(G)
    LG={[]};
    L=0;
end
while ~isempty(G)
    k=k+1;
    nn=G(ceil(rand*size(G,1)),1);
    m1=ismember(G(:,1),nn);
    m2=ismember(G(:,2),nn);
    g=G(m1|m2,:);
    G(m1|m2,:)=[];
    ng=unique(g(:));
    while ~isempty(ng)
        m1=ismember(G(:,1),ng);
        m2=ismember(G(:,2),ng);
        gg=G(m1|m2,:);
        ng=unique(gg(:));
        g=[g;gg];
        G(m1|m2,:)=[];
    end
    LG{k,1}=g;
    L(k,1)=length(unique(g(:)));
    clear g
end