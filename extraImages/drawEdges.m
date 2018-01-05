function im_seg = drawEdges(img,seg)
if ~isa(img,'double')
    img = im2double(img);
end
% [im_r,im_c,~] = size(img);
% im_segr  = zeros(im_r,im_c);
% im_segg  = zeros(im_r,im_c);
% im_segb  = zeros(im_r,im_c);
% img_cr = img(:,:,1);
% img_cg = img(:,:,2);
% img_cb = img(:,:,3);
% for i=1:max(seg(:))
%     idx = seg==i;
%     im_segr(idx) = mean(img_cr(idx));
%     im_segg(idx) = mean(img_cg(idx));
%     im_segb(idx) = mean(img_cb(idx));
% end
% im_seg = cat(3,im_segr,im_segg,im_segb);
im_bound = searchContour(seg);
im_seg = img*0.8+repmat(im_bound,1,1,size(img,3));
im_seg = min(im_seg,1);
end