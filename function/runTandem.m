function allMasks = runTandem(mixture,nChan)

% write a temporal file
dlmwrite('mixture', mixture);

% call the tandem algorithm
cmd = sprintf('!../function/tandem_16000 %d mixture out.tandem', nChan);
eval(cmd);

msk = load(sprintf('out.tandem.%d.mask.dat',nChan));
pc = dlmread(sprintf('out.tandem.%d.pitch.dat',nChan));
pc = pc(:,1:end-1);
nCon = pc(1,1); nFrm = pc(1,2);

% verification
frm1 = 0;
for k=2:nCon+1
    p = pc(k,:);
    frm1 = frm1 + sum(p>0);
end
if frm1 ~= size(msk,1)
    error('Dimensions of pitch and mask do not match!');    
end

% convert to cell structure
count = 0;
msk_st = 1;        
for k=2:nCon+1
    p = pc(k,:);
    if sum(p)
        count = count + 1;
        allMasks{count}.pitch = p;
        allMasks{count}.msk = zeros(128,nFrm);
        st = find(p,1,'first');
        ed = find(p,1,'last');
        dummy = (p(st:ed)==0);                
        msk_ed = msk_st+(ed-st)-sum(dummy);
        idx = st:ed;
        idx = idx(~dummy);
        allMasks{count}.msk(:,idx) = msk(msk_st:msk_ed,:)'>0;
        msk_st = msk_ed+1;
    end
end         

% empirically remove unresonable pitches
tmpMasks = allMasks;
clear allMasks;
mskNum = 0;
for k=1:length(tmpMasks)
    if sum(tmpMasks{k}.pitch>210)<=.3*sum(tmpMasks{k}.pitch>0)
        mskNum = mskNum + 1;
        allMasks{mskNum} = tmpMasks{k};
    end
end

delete out.tandem.128.mask.dat
delete out.tandem.128.pitch.dat