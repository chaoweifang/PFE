function Y = IniGBGaussian(img, paras,nthreads)
if ~exist('nthreads','var') || isempty(nthreads) 
    nthreads = 8;
end
if ~exist('paras','var') || isempty(paras)
    paras.radius1=3;
    paras.spgf1 = 20;
    paras.sigf1 = 10;
    paras.radius2=5;
    paras.spgf2=5;
    paras.sigf2=2;
    paras.scale=1;
    paras.reg = 200;
    paras.lb=100;
    paras.preg=0.8;
end

    %img = imread('../data/12003.jpg');
    %img = im2double(img);
    %tic
    [im_r,im_c,~] = size(img);
    npixels = im_r*im_c;
    img_Lab = applycform(img,makecform('srgb2lab'));
    img_Lab = M_BFilter(img_Lab,paras.radius1,paras.spgf1,paras.sigf1);
    img_rgb = applycform(img_Lab,makecform('lab2srgb'));
    dims = size(img_Lab,3);
    num_center = 2^paras.chs;
    %paras.reg = paras.reg/(321*481) * npixels;
    %paras.lb  = paras.lb/(321*481)  * npixels;
    im_gb_label = mexGBSegMerge(img_rgb*255,1,paras.reg*paras.scale^paras.preg,paras.lb*paras.scale^2,num_center-1);
    assignments = im_gb_label(:);
    num_gb = max(im_gb_label(:));
    initSigmas = zeros(dims,dims,num_gb);
    initMeans = zeros(dims,num_gb);
    X=reshape(img_Lab,npixels,3)';
    N = zeros(num_gb,1);
    %postpdf = zeros(npixels,num_center);
    for i=1:num_gb
      Xk = X(:,assignments==i) ;
      %initWeights(i) = size(Xk,2) / num_center ;
      initSigmas(:,:,i) = cov(Xk')+eye(dims)*10^-3;
      initMeans(:,i) = mean(Xk,2);
      N(i) = size(Xk,2);
    end
    if num_gb ~= num_center-1
        error('Not enough regions from gb')
    end
    %nthreads = 8;
    postpdf = mexGaussianPdf1(X,double(num_gb),nthreads,initMeans,initSigmas,N);
    N=N+eps;
    weights = N/sum(N);
    postpdf = bsxfun(@times, postpdf, weights');
    nc = floor(log(num_center)/log(2));
    k = 2^nc;
    pdfr=zeros(npixels,nc); 
    for i = 1:nc
        len = k/2^i;
        label = (1:2*len:k)';
        label = repmat(label,1,len);
        label = bsxfun(@plus,label,0:len-1);
        label = label(:);
        pdfr(:,i) = sum(postpdf(:,label),2);
    end
    
    Y = zeros(npixels,nc);
    for i=1:nc
        imtemp = reshape(pdfr(:,i),im_r,im_c);
        imtemp = M_BFilter(imtemp,paras.radius2, paras.spgf2, paras.sigf2);
        Y(:,i) = imtemp(:);
    end
    Y = bsxfun(@minus,Y,mean(Y));
end
