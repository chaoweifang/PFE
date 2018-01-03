strimidx = '42049';
teaser_dir='E:\Server5\teaser\60000_100_10';
load(fullfile(teaser_dir,[strimidx '_output.mat']));
vec = output.EigVect;
for i=1:size(vec,2)
    vectemp = vec(:,i);
    im = (reshape(vectemp,size(output.mpb)) - min(vectemp))/(max(vectemp) - min(vectemp));
    imwrite(im,fullfile(teaser_dir,[strimidx '_' num2str(i) '.png']));
end

img = im2double(imread(['E:\Data\BSDS\BSDS500\data\images\val\' strimidx '.jpg']));
%im_contour = im2double(imread(fullfile(teaser_dir,[strimidx '_.mat']))
load(fullfile(teaser_dir,[strimidx '_ucm2.mat']))
figure; imshow(ucm2)
thres = 0.15;
im_contour = ucm2>thres;
im_contour2 = im_contour(2:2:end, 2:2:end);%figure;imshow(im_contour2)
labels2 = bwlabel(ucm2 <= thres);
seg = labels2(2:2:end, 2:2:end);
im_seg = showSegResults(img,seg);
%[gx,gy]=gradient(seg);
%im_bound = gx~=0 + gy~=0;
im_bound = searchContour(seg);
im_seg = im_seg+repmat(im_bound,1,1,size(img,3));
im_seg = min(im_seg,1);%figure;imshow(im_seg)
%im_seg = im_seg+repmat(im_bound,1,1,size(img,3));
%im_seg = min(1,im_seg);
imwrite(im_seg,fullfile(teaser_dir,[strimidx '_seg.png']));
figure;imshow(im_seg)
