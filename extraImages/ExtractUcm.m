addpath('../NcutSplitBregmanCode/benchmarks');
clear

ucmfildir = 'E:\Data\PFE-final-results-20150920\20150920_PFE-MCG-ucm';
if exist(fullfile([ucmfildir '/test_eval'],'eval_cover_img.txt'),'file')
    covering = dlmread(fullfile([ucmfildir '/test_eval'],'eval_cover_img.txt'));
    coveringflag=  true;
else
    coveringflag = false;
end

%ucmfilnames = dir(fullfile(ucmfildir,'*.mat'));
ucmimagedir = [ucmfildir '/ucmimage'];
%imagefilpath = '/home/chiney/Documents/MATLAB/images/test';
%gtfilpath = '/home/chiney/Documents/MATLAB/groundTruth/test';
imagefilpath = 'E:\Data\BSDS\BSDS500\data\images\test';
gtfilpath = 'E:\Data\BSDS\BSDS500\data\groundTruth\test';
contourimagedir = [ucmfildir '/contourimage'];
segimagedir = [ucmfildir '/segimage'];
imgfilnames = dir(fullfile(imagefilpath,'*.jpg'));

if ~exist(ucmimagedir,'dir'), mkdir(ucmimagedir);end
if ~exist(contourimagedir,'dir'), mkdir(contourimagedir);end
if ~exist(segimagedir,'dir'), mkdir(segimagedir);end

for i=1:length(imgfilnames)
    [~,strimidx,~] = fileparts(imgfilnames(i).name);
    ucmfilpath = fullfile(ucmfildir,[strimidx '.mat']);
    load(ucmfilpath);
    imwrite(ucm2,fullfile(ucmimagedir,[strimidx '.png']));
    if coveringflag
        thre = covering(i,2);
    else
        load(fullfile(gtfilpath,[strimidx '.mat']))
        [thresh, cntR, sumR, cntP, sumP, cntR_best,sumRI,sumVOI] = evaluation_reg_ucm_perimage(ucm2, groundTruth, 99);
        Cover = cntR./(sumR + (sumR==0));
        [maxCover,maxIndex] = max(Cover);
        thre = thresh(maxIndex);
    end
    im_contour = ucm2>thre;
    imwrite(im_contour,fullfile(contourimagedir,[strimidx '.png']))

    labels2 = bwlabel(ucm2 <= thre);
    seg = labels2(2:2:end, 2:2:end);
    %figure;imshow(seg/max(seg(:)))
    img = imread(fullfile(imagefilpath,[strimidx '.jpg']));
    img = im2double(img);
    im_seg = showSegResults(img,seg);
    %[gx,gy] = gradient(seg);
    im_bound = searchContour(seg);
    im_seg = im_seg+repmat(im_bound,1,1,size(img,3));
    im_seg = min(im_seg,1);%figure;imshow(im_seg)
    imwrite(im_seg,fullfile(segimagedir,[strimidx '.png']))
end


