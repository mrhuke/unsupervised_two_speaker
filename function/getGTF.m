function f=getGTF(sig,nChan,fs)
% Extract the Gammatone feature (GF)

r = gammatone(sig,nChan,[50 8000],fs);
f = (abs(resample(abs(r)',100,fs))').^(1/3);
