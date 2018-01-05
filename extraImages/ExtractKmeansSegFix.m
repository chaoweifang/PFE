clear;

kmsegfildir = '/home/chiney/Documents/MATLAB/Results/20150428_KmSegBregGMMAffiMPB/seg/breg/fixkm';
segfilnames = dir(fullfile(kmsegfildir,'*.mat'));
imagefilpath = '/home/chiney/Documents/MATLAB/images/test';

scaleimagedir = [kmsegfildir '/scaleimage'];


if ~exist(scaleimagedir,'dir'), mkdir(scaleimagedir);end

for i=1:length(segfilnames)
    [~,strimidx,~] = fileparts(segfilnames(i).name);
    segfilpath = fullfile(kmsegfildir,segfilnames(i).name);
    load(segfilpath);
    
    img = imread(fullfile(imagefilpath,[strimidx '.jpg']));
    img = im2double(img);
    
    %save all scales
    scaleimagedirtemp = [scaleimagedir '/' strimidx];
    mkdir(scaleimagedirtemp);
    for j=1:length(segs)
        seg = double(segs{j});
        im_seg = showSegResults(img,seg);
        [gx,gy] = gradient(seg);
         im_bound = searchContour(seg);
        im_seg = im_seg+repmat(im_bound,1,1,size(img,3));
        im_seg = min(im_seg,1);%figure;imshow(im_seg)
        imwrite(im_seg,fullfile(scaleimagedirtemp,[num2str(j) '.png']))
    end 
end