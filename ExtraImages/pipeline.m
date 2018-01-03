addpath('./benchmarks');
addpath('./lib');
addpath('./Bregman');
strimidx = '101087';
mpbfil = '/home/chiney/Documents/MATLAB/CodeSegmentation/Results/20150418_hidataStep2RGB/mpb/101087.mat';
sgfil = '/home/chiney/Documents/MATLAB/CodeSegmentation/Results/20150418_hidataStep2RGB/sg/101087.mat';
imgdir = '/home/chiney/Documents/MATLAB/CodeSegmentation/data/101087.jpg';
gtFile =strcat('/home/chiney/Documents/MATLAB/CodeSegmentation/dataground/',strimidx,'.mat');
pipelinedir = './pipeline';
if ~exist(pipelinedir,'dir'), mkdir(pipelinedir); end

load(mpbfil);load(sgfil);
img = imread(imgdir);
img = im2double(img);
%img_r = img(:,:,1);
imwrite(img(:,:,1),fullfile(pipelinedir,[strimidx '_rchannel.bmp']))
imwrite(img(:,:,2),fullfile(pipelinedir,[strimidx '_gchannel.bmp']))
imwrite(img(:,:,3),fullfile(pipelinedir,[strimidx '_bchannel.bmp']))

[output,grad] = globalPbBregMK(img,[],'GMM','Step2',sg,mpb);
save(fullfile(pipelinedir,[strimidx '_output.mat']),'output');
save(fullfile(pipelinedir,[strimidx '_grad.mat']),'grad');

%save stage 1 and final embedding
load(fullfile(pipelinedir,'medianresults.mat'));
[im_r,im_c,~] = size(img);
result_stage1 = results(:,:,1);
result_final = results(:,:,end);
for i=1:size(result_stage1,2)
    vectemp = result_stage1(:,i);
    minvectemp = min(vectemp);
    maxvectemp = max(vectemp);
    vectemp = (vectemp-minvectemp)/(maxvectemp-minvectemp);
    imwrite(reshape(vectemp,im_c,im_r)',fullfile(pipelinedir,[strimidx '_stage1_' num2str(i) '.bmp']))
    
    vectemp = result_final(:,i);
    minvectemp = min(vectemp);
    maxvectemp = max(vectemp);
    vectemp = (vectemp-minvectemp)/(maxvectemp-minvectemp);
    imwrite(reshape(vectemp,im_c,im_r)',fullfile(pipelinedir,[strimidx '_final_' num2str(i) '.bmp']))
end

%hierarchical segmentation
ucm2 = contours2ucm(output.gpb_orient, 'doubleSize');
save(fullfile(pipelinedir,[strimidx '_ucm2.mat']),'ucm2');
imwrite(ucm2,fullfile(pipelinedir,[strimidx '_ucm2.bmp']));

load(gtFile);
[thresh, cntR, sumR, cntP, sumP, cntR_best,sumRI,sumVOI] = evaluation_reg_ucm_perimage(ucm2, groundTruth, 99);
Cover = cntR./(sumR + (sumR==0));
[maxCover,maxIndex] = max(Cover);
[maxRI,maxRIIndex] = max(sumRI);
[maxVOI,maxVOIIndex] = min(sumVOI);

fid = fopen(fullfile(pipelinedir,[strimidx '_regionbench.txt']),'w');
score = [thresh(maxIndex) maxCover; thresh(maxRIIndex) maxRI; thresh(maxVOIIndex) maxVOI];
fprintf(fid,'%g %g \n',score');
fclose(fid);

im_contour = ucm2>thresh(maxIndex);
imwrite(im_contour,fullfile(pipelinedir,[strimidx '_contour.bmp']))
figure('name','best covering seg contour');
imshow(ucm2>thresh(maxIndex))
%im_contour2 = im_contour(2:2:end, 2:2:end);figure;imshow(im_contour2)
labels2 = bwlabel(ucm2 <= thresh(maxIndex));
seg = labels2(2:2:end, 2:2:end);
%im_seg = zeros(im_r,im_c,3);
im_seg = showSegResults(img,seg);
%im_b = edge(seg);
%figure;imshow(im_b)
%[bounds,~] = imcontour(seg);
%boundidx = (bounds(2,:)-1)*size(seg,1) + bounds(1,:);
%im_seg2= im_seg;
%im_seg2([boundidx boundidx+numel(seg) boundidx+numel(seg)*2]) = 1.0;
%im_seg = min(im_seg+double(repmat(im_contour(2:2:end,2:2:end),1,1,3)),1);
imwrite(im_seg,fullfile(pipelinedir,[strimidx '_seg.bmp']));
figure('name','best covering seg region');
imshow(im_seg)

%kmeans segmentation
segs_km = calKmeansSegs(output.EigVect,im_r,im_c,5:2:25);
[threshKm, cntRKm, sumRKm, cntPKm, sumPKm, cntR_bestKm,sumRIKm,sumVOIKm] = evaluation_reg_segs_perimage(segs_km,groundTruth);
CoverKm = cntRKm./(sumRKm + (sumRKm==0));
[maxCoverKm,maxIndexKm] = max(CoverKm);
[maxRIKm,maxRIIndexKm] = max(sumRIKm);
[maxVOIKm,maxVOIIndexKm] = min(sumVOIKm);

fid = fopen(fullfile(pipelinedir,[strimidx '_regionbenchKm.txt']),'w');
scoreKm = [threshKm(maxIndexKm) maxCoverKm; threshKm(maxRIIndex) maxRIKm; threshKm(maxVOIIndexKm) maxVOIKm];
fprintf(fid,'%g %g \n',scoreKm');
fclose(fid);

im_seg_Km = showSegResults(img,segs_km{maxIndexKm});
imwrite(im_seg_Km,fullfile(pipelinedir,[strimidx '_seg_kmeans.bmp']));
figure('name','best covering kmeans seg region');
imshow(im_seg_Km)


%save(fullfile(saveucmdir,[strimidx '.mat']),'ucm2');


%evFile4 = fullfile(outDir, strcat(pipelinedir,strimidx, '_ev4.txt'));
% inFile = fullfile(inDir, strcat(pipelinedir,strimidx, '.mat'));
% gtFile = fullfile(gtDir, strcat('/home/chiney/Documents/MATLAB/CodeSegmentation/dataground',strimidx,'.mat'));
% evFile2 = fullfile(outDir, strcat(pipelinedir,strimidx, '_ev2.txt'));
% evFile3 = fullfile(outDir, strcat(pipelinedir,strimidx, '_ev3.txt'));
% 
% evaluation_reg_image(inFile, gtFile, evFile2, evFile3, evFile4, nthresh);
    
    
    