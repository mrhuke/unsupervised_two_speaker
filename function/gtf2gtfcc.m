function [ outData ] = gtf2gtfcc( inData, ccST, ccEND )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

 [chnNum frmNum] = size(inData);
 
 mtx = dctmtx(chnNum);
 
 outData = mtx * inData;
 outData = outData(ccST:ccEND, :);
 
