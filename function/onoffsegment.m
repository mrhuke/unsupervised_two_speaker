function seg = onoffsegment(mixture,params)
% onset/offset segmentation (Hu & Wang'07), input signal will be converted
% to 20 kHz

mixFN = [params.prefix, '.mixture.20k'];    
if params.Srate~=20000     
    mixture20k = resample(mixture,20000,params.Srate);    
end
dlmwrite(mixFN, mixture20k, 'delimiter', ' ');

% cross-channel correlations
% cmd = sprintf('!../function/getCorr %d mixture.20k mixture.corr mixture.evCorr',params.nChan);
% eval(cmd);            
% % remove debug files
% delete test.txt;
% delete test1.txt;

% onset/offset segmentation
corrFN = [params.crossFN, '.', num2str(params.nChan)];
evCorrFN = [params.evCrossFN, '.', num2str(params.nChan)];
segFN = sprintf('%s.%d.seg',params.prefix,params.nChan);
cmd = sprintf('!../function/segmentation %d %s %s %s %s',params.nChan,...
    mixFN,corrFN,evCorrFN,segFN);
eval(cmd);        

seg = binLoad(segFN,128,'int')';
