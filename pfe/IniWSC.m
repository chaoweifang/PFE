function Y = IniWSC(img,paras)
%%input
%I:          t*1        t=pairs of connections
%J:          t*1      
%img:        rs*cs*dims struct chs= number of channels
%chs
%paras       struct {radius2,spgf2,sigf2}

%%output
%Y:             num*chs

    [im_r,im_c,dims]=size(img);
    npixels = im_r*im_c;
    chs = paras.chs;
    %img_Lab = applycform(img,makecform('srgb2lab'));
    [I,J] = cimgnbmap([im_r im_c],5,1);
    X=reshape(img,im_r*im_c,dims);
    %X(:,1) = 0.5*X(:,1);
    Xpdist =sum( (X(I,:)-X(J,:)).^2,2);
    Xmdist = mean(Xpdist)+eps;
    val = exp(-Xpdist/(paras.awsc*Xmdist));
    
    num_center = 2^chs;
    assignments = wspectralclustering(double(I),double(J),val,num_center,num_center);
    num_gb = max(assignments(:));
    initSigmas = zeros(dims,dims,num_gb);
    initMeans = zeros(dims,num_gb);
    %X=reshape(img_Lab,npixels,3)';
    N = zeros(num_gb,1);
    X=X';
    %postpdf = zeros(npixels,num_center);
    for i=1:num_gb
      Xk = X(:,assignments==i) ;
      %initWeights(i) = size(Xk,2) / num_center ;
      initSigmas(:,:,i) = cov(Xk')+eye(dims)*10^-3;
      initMeans(:,i) = mean(Xk,2);
      N(i) = size(Xk,2);
    end
    if num_gb ~= num_center
        error('Not enough regions from gb')
    end
    %nthreads = 8;
    postpdf = mexGaussianPdf1(X,double(num_gb),16,initMeans,initSigmas,N);
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