function vect = im2vect_le(img,edg,param_s,param_ini)

O = edgeOrient(edg,param_s.or1);
e_nms = edgesNmsMex(edg,O,1,5,1.01,4);
img = im2double(img);
img_t=[];
for i=1:size(img,3)
   img_t = cat(3,img_t,img(:,:,i)');
end
lbp = e_nms * param_s.mult_Pb;

nvec = param_ini.chs;
[tx,ty]=size(lbp);
l{1} = zeros(tx + 1, ty);
l{1}(2:end, :) = lbp;
l{2} = zeros(tx, ty + 1);
l{2}(:, 2:end) = lbp;

[val,I,J] = buildW(l{1},l{2}, param_s.dthresh, param_s.ic_gamma);

% if strcmp(param_ini.flag,'g')
%     Yor = IniMtGMM(img_t,param_ini);
% else
%     if strcmp(param_ini.flag,'h')
%         Yor = IniGBGaussian(img_t,param_ini);
%     else if strcmp(param_ini.flag,'s')
%             Yor = IniWSC(img_t,param_ini);
%         else
%             error('Error initialization type!');
%         end
%     end
% end
%figure; ch=3; im = Yor(:,ch); im = (im-min(im))/(max(im)-min(im));im = reshape(im,ty,tx); imshow(im')

[vectmp,~] = NestedBregmanPRW(val,I,J,Yor,param_s.p,param_s.miu,param_s.r,param_s.rex,param_s.rscale,param_s.rlb);
clear val I J Yor;
vect = zeros(tx, ty, nvec);
for v = 1 : nvec,
    temp = reshape(vectmp(:, v), [ty tx])';
    %vect(:,v) = temp(:);
    vect(:, :, v) = temp;%(temp-vmin)/vmm;
end

end