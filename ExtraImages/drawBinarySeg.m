function im_seg = drawBinarySeg(img,seg)
if ~isa(img,'double')
    img = im2double(img);
end
%[im_r,im_c,~] = size(img);

img_cr = img(:,:,1);
img_cg = img(:,:,2);
img_cb = img(:,:,3);
%im_segr  = zeros(im_r,im_c);
%im_segg  = zeros(im_r,im_c);
%im_segb  = zeros(im_r,im_c);

im_bound = searchContour(seg);
im_segr = img_cr;
im_segr(seg==1) = 0.8;
msk_bd = im_bound>0;
im_segr(msk_bd) = 1;
im_segg = img_cg;
im_segg(msk_bd) = 0;
im_segb = img_cb;
im_segb(msk_bd) = 0;
im_seg = cat(3,im_segr,im_segg,im_segb);

%im_seg = im_seg+repmat(im_bound,1,1,size(img,3));
%im_seg = min(im_seg,1);
end