%imidx = 228076
vecfilbase = 'E:\Data\PFE-final-results-20150920\20150922_Km_Seg_BregGMM_AffiNcut_miu10000_r_100_10\vecbreg';
%imgfilbase = 'E:\Data\results\final\VecSegBregGMM\';
imgfilpath = 'E:\Data\BSDS\BSDS500\data\images\test';

imgfilnames = dir(fullfile(vecfilbase,'*.mat'));
 pipelinedir = [vecfilbase '/vecBregGMM'];
if ~exist(pipelinedir,'dir'), mkdir(pipelinedir); end
for i = 1:numel(imgfilnames)
    imidx = imgfilnames(i).name(1:end-4);
   %strimidx=num2str(imidx);
    if exist(fullfile(pipelinedir,[imidx '_vec_' num2str(1) '.png']),'file'), continue; end

    load(fullfile(vecfilbase,[imidx '.mat']));

    img = imread(fullfile(imgfilpath,[imidx '.jpg']));
    [im_c,im_r,~]=size(img);
    %load(fullfile(imagefilpath,[num2str(imidx) '.mat']))
    for j=1:size(vec,2)
        vectemp = vec(:,j);
        minvectemp = min(vectemp);
        maxvectemp = max(vectemp);
        vectemp = (vectemp-minvectemp)/(maxvectemp-minvectemp);
        imwrite(reshape(vectemp,im_c,im_r),fullfile(pipelinedir,[imidx '_vec_' num2str(j) '.png']))
    end
end