function [allMasks params] = runTandem1(mixture, params)

% ensure paralle funning of this function
prefix = [params.workFolder,'/',datestr(now,30),'.tandem'];
while (size(dir([prefix,'*']))>0)
    pause(1);
    prefix = [params.workFolder,'/',datestr(now,30),'.tandem'];
end
params.prefix = prefix;

% write a temporary file
if params.Srate~=20000
    mixture = resample(mixture,20000,params.Srate);
end
params.mixFN = [prefix,'.mixture'];
dlmwrite(params.mixFN, mixture);

% call the tandem algorithm and output response and envelope cross-channel
% correlations
params.tandemFN = prefix;
params.crossFN = [prefix,'.cross'];
params.evCrossFN = [prefix,'.evCross'];
params.engFN = [prefix,'.eng'];
cmd = sprintf('!./tandem %d 20000 .5 %s %s %s %s %s',...
    params.nChan, params.mixFN, params.tandemFN, params.crossFN, params.evCrossFN, params.engFN);
eval(cmd);

msk = load(sprintf('%s.%d.mask.dat',params.tandemFN,params.nChan));
pc = dlmread(sprintf('%s.%d.pitch.dat',params.tandemFN,params.nChan));
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
        allMasks{count}.msk = zeros(params.nChan,nFrm);
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
