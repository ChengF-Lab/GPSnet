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