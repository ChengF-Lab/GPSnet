function Cancer_Module=Cancer_Module_Calculation(Cancer_Type,alpha,cutoff);
%%% Cancer_Module=Cancer_Module_Calculation('LUAD',0.5,0.005);

if nargin < 2
    alpha=0.5;
    cutoff=0.005;
end

%%%%%%% in our paper, we set cutoff=0.5%
str=num2str(100*alpha);
load(['Data_mat/Raw_Module/Raw_Module_',Cancer_Type,'_',str]); %%%%% module

%%%%% top_modules
m=find(Score(:,3)<10);Score(m,:)=[];Module(m)=[]; %% delete the too small module
Score(:,4)=1:length(Score);
Score=sortrows(Score,-2);
score=Score(1:ceil(length(Score)*5/100),:); %%% top 5% module
module_gene=cell2mat(Module(score(:,4))); 
umg=unique(module_gene);

%%%%%%% Cancer_Module
Cancer_Module=[];
for i=1:length(umg)
    l=length(find(module_gene==umg(i)));
    umg(i,2)=l;
end
umg=sortrows(umg,-2);
l=find(umg(:,2)<cutoff*length(score),1,'first');
Cancer_Module=umg(1:l,1);

