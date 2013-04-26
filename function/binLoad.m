function data = binLoad(filename, nCol, type)

fid = fopen(filename,'r');
if (fid>0)
    data = fread(fid,[nCol,inf],type)';
    fclose(fid);
else
    data=[];
end
    
