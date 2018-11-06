function Raw_Module_Generation(Cancer_Type,alpha);

%Cancer_Type='SCLC';
matlabpool('open','local',12);

load(['Data_mat/Cancer_Specific_PPI/',Cancer_Type]);   %%%%% network
load(['Data_mat/Mutation/',Cancer_Type]);  %%%%% mutation
load Data_mat/Gene_Length %%%%% gene length

%%%%%%% eliminate the infulence of the high connected gene
a=ismember(Net,[7316,7273]);
Net(a(:,1)+a(:,2)~=0,:)=[];
[LG,L]=largest_component(Net);
Net=LG{find(L==max(L))};


%%% mutation varify (mutation/gene length)
Mutation=sortrows(Mutation,1);
[a,b]=ismember(Mutation(:,1),Gene_Length(:,1));
Gene_Length(b(a),3)=Mutation(a,2);
Mutation=Gene_Length;
MG=unique(Net(:));
[a,b]=ismember(Mutation(:,1),MG);
MG(b(a),2:3)=Mutation(a,2:3);
PM=MG(:,1);
for i=1:length(MG)
    if MG(i,2)==0
        PM(i,2)=0;
    else
        PM(i,2)=MG(i,3)/MG(i,2);
    end
end

%PM(:,1)=MG(randperm(length(MG)));
[a,G]=ismember(Net,MG(:,1)); %%%%% mutation smoothing
G=sparse([G(:,1);G(:,2)],[G(:,2);G(:,1)],1);
%alpha=0;
F=network_smoothing(PM(:,2),G,alpha,1);
PM(:,2)=F;


parfor II=1:12
    Module_Forming_Process(Net,PM,II,Cancer_Type,alpha);
end

matlabpool close
exit
