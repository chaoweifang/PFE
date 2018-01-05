addpath('./pfe/');
addpath('./external/suitesparse/');
addpath('./external/bsr/')
addpath('./external/vlfeat/')
addpath('./scripts')
addpath('./external/structured_forest/')
addpath('./lib')
addpath('./extraImages/')

clear;clc;

%read image 
img = imread('./tmp/97010.jpg');
%load boundary
load('./tmp/97010_db.mat')

%obtain related parameters
[param_ini,param_s,param_ucm]=get_parameters('s',0.8,'w', 'journal');

%compute initialization, pfe, ucm
[ucm2,~,vecs,~,t]= img2ucm(img,[],[],edg,param_s,param_ucm,param_ini);

[h,w,~]=size(img);
figure('Name','original image')
imshow(img)

figure('Name','pfe')
img_vec=showEmbedding(vecs.vect,w,h);
imshow(img_vec)

figure('Name','segmentation')
seg = (bwlabel(ucm2'<=0.15,4))';
seg = seg(2:2:end,2:2:end);
img_seg=showSegResults(img,seg);
imshow(img_seg)


