function scr= getScr(ind,path,allMasks)

warning off;

dat1 = []; dat2 = [];
for k=1:length(path)
    if path(k)==1
        dat1 = [dat1 allMasks{ind(k)}.dat];
    else if path(k)==2
            dat2 = [dat2 allMasks{ind(k)}.dat];
        end
    end
end
            
grand_me = mean([dat1 dat2],2);
n1 = size(dat1,2);
n2 = size(dat2,2);
Sw = (n1-1)*cov(dat1')+(n2-1)*cov(dat2');
Sb = sum(sum(dat1)~=0)*(mean(dat1,2)-grand_me)*(mean(dat1,2)-grand_me)'...
    +sum(sum(dat2)~=0)*(mean(dat2,2)-grand_me)*(mean(dat2,2)-grand_me)';

if det(Sw)~=0
    scr = trace(Sb*inv(Sw)); 
else
    scr = -1;
end