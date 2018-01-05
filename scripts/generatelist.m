datasetroot='/home2/cwfang/SBD/images';

fils=dir(fullfile(datasetroot,'*.jpg'));

fid=fopen('/home2/cwfang/SBD/test.txt','w');
for i=1:length(fils)
    fprintf(fid,'%s\n',fils(i).name(1:end-4));
end
fclose(fid);


