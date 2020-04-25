
function [PFront PSet]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
PFront=zeros(Xpop,Nobj);
PSet=zeros(Xpop,Nvar);
k=0;

for xpop=1:Xpop
    Dominated=0;
    
    for compare=1:Xpop
        if F(xpop,:)==F(compare,:)
            if xpop > compare
                Dominated=1;
                break;
            end
        else
            if F(xpop,:)>=F(compare,:)
                Dominated=1;
                break;
            end
        end
    end
    
    if Dominated==0
        k=k+1;
        PFront(k,:)=F(xpop,:);
        PSet(k,:)=C(xpop,:);
    end
end
PFront=PFront(1:k,:);
PSet=PSet(1:k,:);
end