function mask = twoSpk_unsupervised(mixture)

% This program implements the algorithm in "An unsupervised approach to cochannel speech separation" by 
% K. Hu and D. L. Wang (submited to IEEE Trans. Audio, Speech, and Lang. Process.). This is an unsupervised
% algorithm for two-speaker separation
%
% mixture: a 16-kHz input two-speaker mixture
% mask   : mask{1} and mask{2} correspond to estimated masks of two speakers

addpath(genpath('../function'));

params = struct('isnr',0,'nChan',128,'Srate',16000,'frate',0.01,'fLen',0.02,...
    'bW',16,'Obj','trace','Pen','linear','p',.1,'Constraint','Soft','lambda',.5,...
    'workFolder','.');

% rescale mixture
mixture = mixture/max(abs(mixture))*2^15;

%% 1. Run a tandem algorithm (Hu & Wang'11) to estimate simultaneous streams
fprintf('Estimating simultaneous streams...\n');
[allMasks params] = runTandem1(mixture,params);
fprintf('Done.\n')

%% 2. Rank simultaneous streams by time
mskNum = length(allMasks);
st = zeros(1,mskNum);
for k=1:mskNum
    st(k) = find(allMasks{k}.pitch>0, 1, 'first');
end
[dummy ind] = sort(st);
outMasks = cell(1,mskNum);
for k=1:mskNum
    outMasks{k} = allMasks{ind(k)};
end
allMasks = outMasks;

%% 3. Extract GFCC features (Shao'07) for each simultaneous stream
gtf = getGTF(mixture,params.nChan,params.Srate); gtf=gtf(:,1:end-1);                
mskNum = length(allMasks);
for k=1:mskNum                            
    msk = allMasks{k}.msk;                        
    mgtf = gtf.*msk;
    CCData = gtf2gtfcc(mgtf, 1, 30);
    allMasks{k}.dat = CCData(:,allMasks{k}.pitch>0);
end

%% 4. Perform beam search to group voiced simultaneous streams
fprintf('\nGrouping voiced speech...\n');
estLab = BeamSearch(allMasks, params);
fprintf('Done.\n')

% produce voiced masks
nFrame = size(allMasks{1}.msk,2);
vMask{1}=zeros(params.nChan,nFrame);
vMask{2}=zeros(params.nChan,nFrame);
pc{1}=[]; pc{2}=[];
for i=1:length(estLab)
    if estLab(i)
        vMask{1} = sign(vMask{1} + allMasks{i}.msk);
        pc{1} = [pc{1}; allMasks{i}.pitch];
    else
        vMask{2} = sign(vMask{2} + allMasks{i}.msk);
        pc{2} = [pc{2}; allMasks{i}.pitch];
    end
end

%% 5. Onset/offset segmentation to generate unvoiced segments
fprintf('\nExtracting unvoiced speech segments...\n');
ooSegs = onoffsegment(mixture,params);
uvSegs = ooSegs;
uvSegs(:,(sum(pc{1},1)>0 & sum(pc{2},1)>0)) = 0;
uvSegs(:,(sum(pc{1},1)==0 & sum(pc{2},1)==0)) = 0;
uvSegs = remvOverlap(uvSegs,vMask{1},vMask{2});  
fprintf('Done.\n')

%% 6. Group unvoiced-voiced (UV) and unvoiced-unvoiced (UU) segments
fprintf('\nGrouping unvoiced speech...');
% group unvoiced speech in UV portions
[nChan,nFrame] = size(vMask{1});
allOneMsk1 = ones(nChan,nFrame); allOneMsk1(:,sum(pc{1},1)==0)=0;
allOneMsk2 = ones(nChan,nFrame); allOneMsk2(:,sum(pc{2},1)==0)=0;
cMsk1 = (1-vMask{1}).*allOneMsk1.*(1-allOneMsk2);
cMsk2 = (1-vMask{2}).*allOneMsk2.*(1-allOneMsk1);
onoffMsk1 = zeros(nChan,nFrame);
onoffMsk2 = zeros(nChan,nFrame);
g = gammatone(mixture, params.nChan, [50,8000], params.Srate);
eng = cochleagram(g, params.fLen*params.Srate);
for k=1:max(max(uvSegs))
    seg = uvSegs==k;
    overlap1 = sum(sum(seg.*cMsk1.*eng));
    overlap2 = sum(sum(seg.*cMsk2.*eng));
    if overlap1>0 || overlap2>0
        if  overlap2>=overlap1
            onoffMsk1 = sign(onoffMsk1 + seg);
        else
            onoffMsk2 = sign(onoffMsk2 + seg);
        end
    end
end
uvMask{1} = sign(onoffMsk1);
uvMask{2} = sign(onoffMsk2);
uvMask{1}(:,~(sum(pc{1},1)==0 & sum(pc{2},1)>0))=0;
uvMask{2}(:,~(sum(pc{2},1)==0 & sum(pc{1},1)>0))=0;          

% group unvoiced speech in UU portions                                
uuSegs = ooSegs;
uuSegs(:,(sum(pc{1},1)>0 | sum(pc{2},1)>0)) = 0;
uuSegs = remvOverlap(uuSegs, vMask{1}, vMask{2});    
uuMask{1} = .5*sign(uuSegs);
uuMask{2} = uuMask{1};

fprintf('Done.\n')

%% produce overall masks and separated signals
mask{1} = vMask{1}+uvMask{1}+uuMask{1};
mask{2} = vMask{2}+uvMask{2}+uuMask{2};
sig{1} = synthesis(mixture,mask{1},[50,8000],params.fLen*params.Srate,params.Srate);
sig{2} = synthesis(mixture,mask{2},[50,8000],params.fLen*params.Srate,params.Srate);
wavwrite(sig{1}/max(abs(sig{1})), params.Srate, 'segregated.1.wav');
wavwrite(sig{2}/max(abs(sig{2})), params.Srate, 'segregated.2.wav');
wavwrite(mixture/max(abs(mixture)), params.Srate, 'mixture.wav');

% remove temporary files
files = dir([params.prefix,'*']);
folder = fileparts(params.prefix);
for k=1:length(files)
    delete([folder,'/',files(k).name]);
end
