function allSegs = remvOverlap(segments,mask1,mask2)

% removing overlappings
allSegs = zeros(size(segments));
nSeg = max(max(segments));
nNewSeg = 0;
for n = 1:nSeg
    seg_n = (segments == n);            
    idx1 = find( (mask1 & seg_n)> 0);
    idx2 = find( (mask2 & seg_n)> 0);
    seg_n(idx1) = 0;
    seg_n(idx2) = 0;
    %find connected regions
    seg_n = bwlabel(seg_n, 4);
    nSubSeg = max(max(seg_n));
    for i = 1:nSubSeg
        newSeg = (seg_n == i);
        tmpSeg.msk = newSeg;
        tmpSeg.pitch = sum(newSeg);
        nNewSeg = nNewSeg + 1;                        
        allSegs = allSegs + nNewSeg*newSeg;
    end
end        