function Module_Treat;


for alpha=0:0.05:1
str=num2str(100*alpha);

Cancer_Type='THCA';

module=[];
score=[];
for i=1:12
    load(['Raw_Module_Score/RandomSeed_Connectivity_Score_Module_',Cancer_Type,'_',str,'_',num2str(i)])
    if exist('Module')
        module=[module;Module];
        score=[score;Score];
        clear Module Score
    end
end
m=find(score(:,1)==0);
score(m,:)=[];
module(m,:)=[];
Module=module;
Score=score;
disp(length(Score))
save(['Data_mat\Raw_Module\','PPI_',Cancer_Type,'_',str],'Module','Score');

end
