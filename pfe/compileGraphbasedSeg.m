%mex -outdir ../lib mexGraphbasedSegFix.cpp

% mex -outdir ../lib mexGraphbasedSegFix.cpp CXXFLAGS='$CXXFLAGS -std=c++0x'
% mex -outdir ../lib mexGraphbasedSegFix1.cpp CXXFLAGS='$CXXFLAGS -std=c++0x'
% mex -outdir ../lib mexGraphbasedSegFix2.cpp CXXFLAGS='$CXXFLAGS -std=c++0x'
 mex mexGBSegMerge.cpp CXXFLAGS='$CXXFLAGS -std=c++0x'
%mex mexGBSegMerge.cpp
%img =im2double( imread('E:\Data\BSDS\BSDS500\data\images\val\86016.jpg'));
%im_seg = mexGBSegMerge(img,1,0.1,100,16);
%figure;imshow(label2rgb(im_seg))
% addpath('../ExtraImages/');
% addpath('../filtering/')
% addpath('../vl/')
% addpath('../mtgmm/')
%addpath('E:\toolbox\gco\matlab\bin')
%addpath('E:\toolbox\gco\matlab')
% 
% img_lab = applycform(img,makecform('srgb2lab'));
% img_bifil = M_BFilter(img_lab,5,20,10);
% img_rgb = applycform(img_bifil,makecform('lab2srgb'));
% num_center = 16;
% im_gb_label = mexGraphbasedSegFix(img_rgb*255,1,100,200,16);
% figure;
% im = double(im_gb_label);
% imshow(im/max(im(:)))
% X=reshape(img_lab,321*481,3)';
% 
% assignments = im_gb_label(:);
% num_gb = max(im_gb_label(:));
% initSigmas = zeros(size(X,1),size(X,1),num_center);
% initMeans = zeros(size(X,1),num_center);
% dims = size(X,1);
% for i=1:num_gb
%   Xk = X(:,assignments==i) ;
%   %initWeights(i) = size(Xk,2) / num_center ;
%   initSigmas(:,:,i) = cov(Xk')+eye(dims)*10^-3;
%   initMeans(:,i) = mean(Xk,2);
% end
% if num_gb < num_center
%     [initMeans1, assignments1] = vl_kmeans(X,num_center-num_gb,'Initialization', 'plusplus','MaxNumIterations',500,'NumRepetitions',10);
%     for i=num_gb+1:num_center
%         Xk = X(:,assignments1==i-num_gb) ;
%         initMeans(:,i) = initMeans1(:,i-num_gb);
%         initSigmas(:,:,i) = cov(Xk')+eye(dims)*10^-3;
%     end
% end
% [clusters,~,phi,postpdf,~] = mexMtGMM(X,num_center,8,initMeans,initSigmas);
% [~,assignGMM] = max(postpdf,[],2);
% im_gmm = reshape(double(assignGMM),321,481);
% label_gmm=denoisingAreas(im_gmm,clusters,3);
% im_gmm_seg = drawEdges(img,im_gmm);
% figure
% imshow(im_gmm_seg)
% 
% nlabels = size(postpdf,2);
% [height,width,~]=size(img);
% W = neighbormap(img_lab,2,30)*20;
% npixels = height*width;
% h_gco = GCO_Create(npixels,nlabels);
% GCO_SetDataCost(h_gco,-log(postpdf+eps)');
% GCO_SetNeighbors(h_gco,W);
% GCO_Swap(h_gco,200);
% label_results = GCO_GetLabeling(h_gco);
% GCO_Delete(h_gco);
% gc_label_img=reshape(double(label_results),321,481);
% gc_seg_img = drawEdges(img,gc_label_img);
% 
% figure
% imshow(gc_seg_img)
% 
% 
% 
% % [initMeans, assignments] = vl_kmeans(X,16,'Initialization', 'plusplus','MaxNumIterations',1000,'NumRepetitions',10);
% %  imidx = reshape(assignments,321,481);
% %        figure;
% %        imshow(double(imidx)/16)
% %        colormap(hsv)
% % figure;
% % imshow(double(im_gb_label)/58)
% % colormap(hsv)
% 
% % addpath('../ExtraImages/')
% %im_label_gbd=denoisingAreas(im_gb_label,clusters,3);
% im_gb_seg = showSegResults(img,double(im_gb_label));
% 
% 
% figure
% imshow(im_gb_seg)