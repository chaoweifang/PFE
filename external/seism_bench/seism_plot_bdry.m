function seism_plot_reg(inDir)
eval_dir = fullfile(inDir,'eval_region.txt');
if ~exist(eval_dir,'file'), error('wrong input dir'); end

fid = fopen(eval_dir);
str_fop = fgetl(fid);
fop_line = fgetl(fid);
fop = cell2mat(textscan(fop_line,'%f %f %f'));
str_sc = fgetl(fid);
sc_line = fgetl(fid);
sc = cell2mat(textscan(sc_line,'%f %f %f'));
str_pri = fgetl(fid);
pri_line=fgetl(fid);
pri = cell2mat(textscan(pri_line,'%f %f %f'));
str_voi = fgetl(fid);
voi_line = fgetl(fid);
voi = cell2mat(textscan(voi_line,'%f %f %f'));

disp([str_fop ' [' num2str(fop(1)) ']' 'ODS=' num2str(fop(2)) ' OIS=' num2str(fop(3))])
disp([str_sc ' [' num2str(sc(1)) ']' 'ODS=' num2str(sc(2)) ' OIS=' num2str(sc(3))])
disp([str_pri ' [' num2str(pri(1)) ']' 'ODS=' num2str(pri(2)) ' OIS=' num2str(pri(3))])
disp([str_voi ' [' num2str(voi(1)) ']' 'ODS=' num2str(voi(2)) ' OIS=' num2str(voi(3))])

fclose(fid);
%evals = fread(fid);
end