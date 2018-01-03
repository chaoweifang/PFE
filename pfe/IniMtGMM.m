function Y = IniMtGMM(img,paras,nthreads)
if nargin < 3
    nthreads = 16;
end
if nargin < 2
    paras.radius1=5;
    paras.spgf1 = 20;
    paras.sigf1 = 10;
    paras.radius2=5;
    paras.spgf2=5;
    paras.sigf2=2;
end
    img_Lab = applycform(img,makecform('srgb2lab'));
    img_Lab = M_BFilter(img_Lab,paras.radius1,paras.spgf1,paras.sigf1);
    img_rgb = applycform(img_Lab,makecform('lab2srgb'));
    [im_r,im_c,~] = size(img);
    Y0_rgb = reshape(img_rgb,im_r*im_c,3);
    num_center = 2^paras.chs;
    X = Y0_rgb';
    [initMeans, assignments] = vl_kmeans(X,num_center,'Initialization', 'plusplus','MaxNumIterations',100,'NumRepetitions',3);
    initSigmas = zeros(size(X,1),size(X,1),num_center);
    initWeights = zeros(1,num_center);
    for i=1:num_center
      Xk = X(:,assignments==i) ;
      initWeights(i) = size(Xk,2) / num_center ;
      if size(Xk,1) == 0 || size(Xk,2) == 0
        initSigmas(:,:,i) = cov(X')+eye(3)*10^-3;
      else
        initSigmas(:,:,i) = cov(Xk')+eye(3)*10^-3;
      end
    end
    [~,~,phi,postpdf,~] = mexMtGMM(X,16,nthreads,initMeans,initSigmas);
    m = im_r*im_c;
    nc = floor(log(num_center)/log(2));
    k = 2^nc;
    %phi = GMModel.PComponents;
    [~,idxsort]  = sort(phi);
    pdf_wr = postpdf(:,idxsort);
    pdfr=zeros(m,nc);
    for i = 1:nc
        len = k/2^i;
        label = (1:2*len:k)';
        label = repmat(label,1,len);
        label = bsxfun(@plus,label,0:len-1);
        label = label(:);
        pdfr(:,i) = sum(pdf_wr(:,label),2);
    end
    
    Y = zeros(m,nc);
    for i=1:nc
        imtemp = reshape(pdfr(:,i),im_r,im_c);
        imtemp = M_BFilter(imtemp,paras.radius2,paras.spgf2,paras.sigf2);
        Y(:,i) = imtemp(:);
    end
    
    Y = bsxfun(@minus,Y,mean(Y));
end
