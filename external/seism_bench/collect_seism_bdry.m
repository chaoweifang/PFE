function collect_seism_bdry(outDir,thresh)

fname = fullfile(outDir, 'eval_bdry.txt');
if exist(fname,'file')
   return 
end

S = dir(fullfile(outDir,'*_ev.txt'));
fb_value=[];
fb_fprec=[];
fb_frec = [];
%voi = [];
%%pri = [];
%sc = [];
for i=1:numel(S)
    evFile = fullfile(outDir,S(i).name);
    vals = dlmread(evFile);
    fb_value = cat(2,fb_value,vals(:,1));
    fb_fprec = cat(2,fb_fprec,vals(:,2));
    fb_frec = cat(2,fb_frec,vals(:,3));
    %if i==1
    %    thresh = vals(:,1);
    %end
    %pri = cat(2,pri,vals(:,4));
    %voi = cat(2,voi,vals(:,5));
    %sc = cat(2,sc,vals(:,6));
end
fb_fprec_ods = mean(fb_fprec,2);
fb_frec_ods = mean(fb_frec,2);
fb_value_ods = 2*fb_fprec_ods.*fb_frec_ods./(fb_fprec_ods+fb_frec_ods);
[fb_ods,fb_ods_tidx]=max(fb_value_ods,[],1);

[Ru, indR] = unique(fb_frec_ods);
Pu = fb_fprec_ods(indR);
Ri = 0:0.01:1;
if numel(Ru)>1
     P_int1 = interp1(Ru, Pu, Ri);
     P_int1(isnan(P_int1)) = 0;
     Area_PR = sum(P_int1)*0.01;
else
     Area_PR = 0;
     P_int1 = zeros(1,numel(Ri));
end

[fb_value_img,fb_value_tidx] = max(fb_value,[],1);
idx = sub2ind(size(fb_value),fb_value_tidx,1:numel(S));
fb_fprec_img = mean(fb_fprec(idx));
fb_frec_img  = mean(fb_frec(idx));
fb_ois = 2*fb_fprec_img*fb_frec_img/(fb_fprec_img+fb_frec_img);

% sc_thr = mean(sc,2);
% [sc_ods,sc_ods_tidx] = max(sc_thr);
% [sc_img,sc_img_tidx] = max(sc,[],1);
% sc_ois = mean(sc_img);
% %disp(sc_ods)
% %disp(sc_ois)
% 
% voi_thr = mean(voi,2);
% [voi_ods,voi_ods_tidx] = min(voi_thr);
% [voi_img,voi_img_tidx] = min(voi,[],1);
% voi_ois = mean(voi_img);
% %disp(pri_ods)
% %disp(pri_ois)
% 
% 
% pri_thr = mean(pri,2);
% [pri_ods,pri_ods_tidx] = max(pri_thr);
% [pri_img,pri_img_tidx] = max(pri,[],1);
% pri_ois = mean(pri_img);

fname=fullfile(outDir,'eval_bdry.txt');
fid = fopen(fname,'w');
fprintf(fid,'%s:\n','Fb');
fprintf(fid,'%10g %10g %10g %10g\n', thresh(fb_ods_tidx), fb_ods, fb_ois, Area_PR);
% fprintf(fid,'%s:\n','SC');
% fprintf(fid,'%10g %10g %10g\n', thresh(sc_ods_tidx), sc_ods, sc_ois);
% fprintf(fid,'%s:\n','PRI');
% fprintf(fid,'%10g %10g %10g\n', thresh(pri_ods_tidx), pri_ods, pri_ois);
% fprintf(fid,'%s:\n','VOI');
% fprintf(fid,'%10g %10g %10g\n', thresh(voi_ods_tidx), voi_ods, voi_ois);
fclose(fid);
% 
% fname=fullfile(outDir,'eval_voi_th.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10g %10g\n',[thresh voi_thr]');
% fclose(fid);
% fname=fullfile(outDir,'eval_voi_img.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(voi_img_tidx) voi_img']');
% fclose(fid);
% 
% fname=fullfile(outDir,'eval_pri_th.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10g %10g\n',[thresh pri_thr]');
% fclose(fid);
% fname=fullfile(outDir,'eval_pri_img.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(pri_img_tidx) pri_img']');
% fclose(fid);
% 
% fname=fullfile(outDir,'eval_cover_th.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10g %10g\n',[thresh sc_thr]');
% fclose(fid);
% fname=fullfile(outDir,'eval_cover_img.txt');
% fid = fopen(fname,'w');
% fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(sc_img_tidx) sc_img']');
% fclose(fid);

fname=fullfile(outDir,'eval_fb_th.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[thresh fb_value_ods]');
fclose(fid);
fname=fullfile(outDir,'eval_fb_img.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10d %10g %10g\n',[(1:numel(S))' thresh(fb_value_tidx) fb_value_img']');
fclose(fid);
fname=fullfile(outDir,'eval_fb_pr.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g %10g\n',[thresh fb_fprec_ods fb_frec_ods]');
fclose(fid);
fname=fullfile(outDir,'eval_fb_prcurve.txt');
fid = fopen(fname,'w');
fprintf(fid,'%10g %10g\n',[Ri;P_int1]);
fclose(fid);

%disp(voi_ods)
%disp(voi_ois)


end