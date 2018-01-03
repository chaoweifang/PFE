function Y = IniGBC(data,pair_in_neigh,paras)
%%input
%data:          num*dim
%pair_in_neigh: t*2      t=pairs of connections
%paras:         reg= regularization parameter for region shape
%               ms = minimum region size
%               chs= number of channels
%%output
%Y:             num*chs

idx1 = uint32(pair_in_neigh(:,1));
idx2 = uint32(pair_in_neigh(:,2));

[~,postpdf] = mexGraphbasedClusterFix(idx1,idx2,single(data)',paras.reg,paras.ms,2^paras.chs);
postpdf = postpdf';
num_center = size(postpdf,2);
nc = floor(log(num_center)/log(2));
k = 2^nc;
ndata = size(data,1);
pdfr=zeros(ndata,nc); 
for i = 1:nc
    len = k/2^i;
    label = (1:2*len:k)';
    label = repmat(label,1,len);
    label = bsxfun(@plus,label,0:len-1);
    label = label(:);
    pdfr(:,i) = sum(postpdf(:,label),2);
end

%Y = zeros(npixels,nc);
%for i=1:nc
%    imtemp = reshape(pdfr(:,i),im_r,im_c);
%    imtemp = M_BFilter(imtemp,paras.radius2, paras.spgf2, paras.sigf2);
%    Y(:,i) = imtemp(:);
%end

Y = bsxfun(@minus,pdfr,mean(pdfr));
    
end