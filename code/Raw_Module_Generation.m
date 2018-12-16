function Raw_Module_Generation(Cancer_Type,alpha);

%%%% Raw_Module_Generation('LUAD',0.5);

if nargin < 2
    alpha=0.5;
end

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
F=network_smoothing(PM(:,2),G,alpha,1);
PM(:,2)=F;

%%%% generate raw modules (the number of raw module is LL) 
LL=60000;
Module_Forming_Process(Net,PM,Cancer_Type,alpha,LL);


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


%%%%%% network smoothing using the network spreading process
function F=network_smoothing(Mutation,G,alpha,flag);
 
Mutation=Mutation/sum(Mutation)*length(Mutation);
F=Mutation;
Y=Mutation;
sn=sum(G,2);
D=sparse(diag(1./sn));
W1=G*D;   %%%% md
delta=100;
k=0;
x=F;
switch flag
    case 1  %%%%% diffusion
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*W1*F0+alpha*Y;
            delta=sum((F-F0).^2);
        end  
    case 2  %%%% conduction
        F=F';
        Y=Y';
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*F0*W1+alpha*Y;
            delta=sum((F-F0).^2);
        end
        F=F';
    case 3
        W1=G;
        [a,b]=find(W1~=0);
        for i=1:length(a)
            W1(a(i),b(i))=1/sqrt(sn(a(i))*sn(b(i)));
        end
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*W1*F0+alpha*Y;
            delta=sum((F-F0).^2);
        end 
end
F=F*sum(Mutation)/length(Mutation);


%%%%%%%% generate the raw module
function Module_Forming_Process(Net,PM,Cancer_Type,alpha,LL);
%%%%%% module merging, zscore=sum(pm-<pm>)/sqrt(k)
%%%%%% seclect probability pm*fi(1-p)/sum(pm*fi(1-p))
tic;
Str_alpha=num2str(100*alpha);

Average_S=mean(PM(:,2));
Gene_List=unique(Net(:));
N=length(Gene_List);
[a,b]=ismember(Net,Gene_List);
Net=sparse([b(:,1);b(:,2)],[b(:,2),b(:,1)],1);
Degree=sum(Net,2);
Module{LL,1}=[];
Score=zeros(LL,3);
Seed_Gene=1:length(Gene_List);
for i=1:LL
    seed=Seed_Gene(ceil(rand*length(Seed_Gene)));
    node=seed;
    score=sum(PM(node,2)-Average_S)/sqrt(1);
    for j=2:1000
        [a,b]=find(Net(:,node)==1);
        node_extend=setdiff(a,node);
        p_extend=[];
        score1=(sum(PM(node,2)-Average_S)+PM(node_extend,2)-Average_S)/sqrt(j);
        m=find(score1>1.01*score);
        if isempty(m)
            break
        else
            node_extend=node_extend(m);
            score1=score1(m);
            for k=1:length(node_extend)
                ks_extend=length(intersect(find(Net(node_extend(k),:)==1),node));
                k_extend=Degree(node_extend(k));
                p_extend(k,1)= sum(hygepdf(ks_extend:k_extend,N,j-1,k_extend));
            end
            node_p_s=[node_extend p_extend score1];            
            node_p_s(find(node_p_s(:,2)>0.05),:)=[];
            if isempty(node_p_s)
                break
            else    
                node_p_s=sortrows(node_p_s,-3);
                x=PM(node_p_s(:,1),2)/sum(PM(node_p_s(:,1),2));
                x=[0;cumsum(x)];
                xx=find(x<rand,1,'last');
                node=[node;node_p_s(xx,1)];
                score=node_p_s(xx,3);
            end
        end
    end
    Module{i,1}=Gene_List(node);
    Score(i,1:3)=[Gene_List(seed) score length(node)];
    if mod(i,1)==0
        toc;
        disp(i)
        save(['Data_mat/Raw_Module/Raw_Module_',Cancer_Type,'_',Str_alpha,'.mat'],'Module','Score')
    end
end