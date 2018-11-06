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
        