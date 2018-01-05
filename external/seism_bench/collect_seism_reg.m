function collect_seism_reg(outDir,thresh)

fname = fullfile(outDir, 'eval_fop_pr.txt');
if exist(fname,'file')
   return 
end

S = dir(fullfile(outDir,'*_ev.txt'));
fop_value=[];
fop_fprec=[];
fop_frec = [];
voi = [];
pri = [];
sc = [];
for i=1:numel(S)
    
    evFile = fullfile(outDir,S(i).name);
    vals = dlmread(evFile);
    fop_value = cat(2,fop_value,vals(:,1));
    fop_fprec = cat(2,fop_fprec,vals(:,2));
    fop_frec = cat(2,fop_frec,vals(:,3));
    pri = cat(2,pri,vals(:,4));
    voi = cat(2,voi,vals(:,5));
    sc = cat(2,sc,vals(:,6));
end
fop_fprec_ods = mean(fop_fprec,2);
fop_frec_ods = mean(fop_frec,2);
fop_value_ods = 2*fop_fprec_ods.*fop_frec_ods./(fop_fprec_ods+fop_frec_ods);
[fop_ods,fop_ods_tidx]=max(fop_value_ods,[],1);

[fop_value_img,fop_value_tidx] = max(fop_value,[],1);
idx = sub2ind(size(fop_value),fop_value_tidx,1:numel(S));
fop_fprec_img = mean(fop_fprec(idx));
fop_frec_img  = mean(fop_frec(idx));
fop_ois = 2*fop_fprec_img*fop_frec_img/(fop_fprec_img+fop_frec_img);

sc_thr = mean(sc,2);
[sc_ods,sc_ods_tidx] = max(sc_thr);
[sc_img,sc_img_tidx] = max(sc,[],1);
sc_ois = mean(sc_img);
%disp(sc_ods)
%disp(sc_ois)

voi_thr = mean(voi,2);
[voi_ods,voi_ods_tidx] = min(voi_thr);
[voi_img,voi_img_tidx] = min(voi,[],1);
voi_ois = mean(voi_img);
%disp(pri_ods)
%disp(pri_ois)


pri_thr = mean(pri,2);
[pri_ods,pri_ods_tidx] = max(pri_thr);
[pri_img,pri_img_tidx] = max(pri,[],1);
pri_ois = mean(pri_img);

fname=fullfile(outDir,'eval_region.txt');
fid = fopen(fname,'w');
fprintf(fid,'%s:\n','FOP');
fprintf(fid,'%10g %10g %10g\n', thresh(fop_ods_tidx), fop_ods, fop_ois);
fprintf(fid,'%s:\n','SC');
fprintf(fid,'%10g %10g %10g\n', thresh(sc_ods_tidx), sc_ods, sc_ois);
fprintf(fid,'%s:\n','PRI');
fprintf(fid,'%10g %10g %10g\n', thresh(pri_ods_tidx), pri_ods, pri_ois);
fprintf(fid,'%s:\n','VOI');
fprintf(fid,'%10g %10g %10g\n', thresh(voi_ods_tidx), voi_ods, voi_ois);
fclose(fid);

fname=fullfile(outDir,'eval_voi_th.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[thresh voi_thr]');
fclose(fid);
fname=fullfile(outDir,'eval_voi_img.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(voi_img_tidx) voi_img']');
fclose(fid);

fname=fullfile(outDir,'eval_pri_th.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[thresh pri_thr]');
fclose(fid);
fname=fullfile(outDir,'eval_pri_img.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(pri_img_tidx) pri_img']');
fclose(fid);

fname=fullfile(outDir,'eval_cover_th.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[thresh sc_thr]');
fclose(fid);
fname=fullfile(outDir,'eval_cover_img.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(sc_img_tidx) sc_img']');
fclose(fid);

fname=fullfile(outDir,'eval_fop_th.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[thresh fop_value_ods]');
fclose(fid);
fname=fullfile(outDir,'eval_fop_img.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(fop_value_tidx) fop_value_img']');
fclose(fid);
fname=fullfile(outDir,'eval_fop_pr.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g %10g\n',[thresh fop_fprec_ods fop_frec_ods]');
fclose(fid);

%disp(voi_ods)
%disp(voi_ois)


end