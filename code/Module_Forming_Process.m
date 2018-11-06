function Module_Forming_Process(Net,PM,II,Cancer_Type,alpha);

%%%%%% module merging, zscore=sum(pm-<pm>)/sqrt(k)
%%%%%% seclect probability pm*fi(1-p)/sum(pm*fi(1-p))

tic;

Str_alpha=num2str(100*alpha);

ctime=datestr(now,30);
tseed=str2num(ctime((end-5):end));
rand('seed',tseed*II);

Average_S=mean(PM(:,2));
Gene_List=unique(Net(:));
N=length(Gene_List);
[a,b]=ismember(Net,Gene_List);
Net=sparse([b(:,1);b(:,2)],[b(:,2),b(:,1)],1);
Degree=sum(Net,2);
LL=5000;
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
    if mod(i,100)==0
        toc;
        disp(i)
        save(['Raw_Module_Score/RandomSeed_Connectivity_Score_Module_',Cancer_Type,'_',Str_alpha,'_',num2str(II),'.mat'],'Module','Score')
    end
end
