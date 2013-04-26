function estLab = BeamSearch(allMasks, params)

mskNum = length(allMasks);

%% Start from the two simultaneous streams with most overlapping in time
allPitch = [];
for k=1:mskNum
    allPitch = [allPitch; sign(allMasks{k}.pitch)];
end
overlap = allPitch*allPitch';
overlap = overlap.*(1-eye(length(allMasks)));
if any(any(overlap))
    [tmp1,tmp2] = find(overlap==max(max(overlap)));                    
    i=tmp1(1); j=tmp2(1);    
else % no overlapping
    i=1;j=2;
end                                           

%% Beam search
estLab = zeros(1,mskNum);
estLab(i)=1; estLab(j)=0;
% exclude the previous two simultaneous streams
clear newMasks;
ct = 1; ind=zeros(1,mskNum-2);
for k = 1:mskNum
    if k~=i && k~=j
        newMasks{ct} = allMasks{k};
        ind(ct) = k;
        ct = ct + 1;
    end
end
% search
grpMasks{1}=allMasks{i}; grpMasks{2}=allMasks{j};    
allLab = [1,0];
for k=1:length(newMasks)
    cW = size(allLab,1);
    % hypothesize the kth simultaneous stream
    grpMasks{end+1} = newMasks{k};
    allLab = [allLab zeros(cW,1); allLab ones(cW,1)];    
    % pruning
    if size(allLab,1)>params.bW || k==length(newMasks)
        scr = objfun(params.Obj, params.Pen, allLab, grpMasks, params.p, params.lambda);
        [sScr sInd]= sort(scr,'descend');
        curInd = sInd( 1:min(params.bW,length(sInd)) );
        allLab = allLab(curInd,:);
    end
end

%% organize labels
estLab = zeros(1,size(allLab,2));
estLab(i) = 1;
estLab(j) = 0;
tmpInd = 1:size(allLab,2);
estLab(tmpInd~=i & tmpInd~=j) = allLab(1,3:end);