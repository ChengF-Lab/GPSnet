function Module_Best=Module_Selected(Cancer_Type,alpha,cutoff);


%%%%%%% in our paper, we set cutoff=0.5%
str=num2str(100*alpha);

%load(['Data_mat/Cancer_Specific_PPI/',Cancer_Type]);   %%%%% network
%load(['Data_mat/Mutation/',Cancer_Type]);  %%%%% mutation
load(['Data_mat/Raw_Module/',Cancer_Type,'_',str]); %%%%% module

%%%%% 
a=ismember(Net,[7316,7273]);
Net(a(:,1)+a(:,2)~=0,:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};

%%%%% best module
m=find(Score(:,3)<10);Score(m,:)=[];Module(m)=[]; %% delete the too small module
module=Module;
Score(:,4)=1:length(Score);
Score=sortrows(Score,-2);
score=Score(1:ceil(length(Score)*5/100),:); %%% top 1% module
module_gene=cell2mat(Module(score(:,4))); 
umg=unique(module_gene);

%%%%%%% using cutoff
Module_Best=[];
for i=1:length(umg)
    l=length(find(module_gene==umg(i)));
    umg(i,2)=l;
end
umg=sortrows(umg,-2);
l=find(umg(:,2)<cutoff*length(score),1,'first');
Module_Best=umg(1:l,1);

