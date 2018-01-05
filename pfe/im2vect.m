function vect = im2vect(Yor,edg,param_s)

O = edgeOrient(edg,param_s.or1);
e_nms = edgesNmsMex(edg,O,1,5,1.01,4);

lbp = e_nms * param_s.mult_Pb;

nvec = size(Yor,2);
[tx,ty]=size(lbp);
l{1} = zeros(tx + 1, ty);
l{1}(2:end, :) = lbp;
l{2} = zeros(tx, ty + 1);
l{2}(:, 2:end) = lbp;

[val,I,J] = buildW(l{1},l{2}, param_s.dthresh, param_s.ic_gamma);


[vectmp,~] = NestedBregmanPRW(val,I,J,Yor,param_s.p,param_s.miu,param_s.r,param_s.rex,param_s.rscale,param_s.rlb);
clear val I J Yor;
vect = zeros(tx, ty, nvec);
%vectre = zeros(tx,ty,nvec);
for v = 1 : nvec,
    temp = reshape(vectmp(:, v), [ty tx])';
    %vect(:,v) = temp(:);
    vect(:, :, v) = temp;%(temp-vmin)/vmm;
    %vmin = min(temp(:)); vmax=max(temp(:)); vmm = vmax-vmin;
    %tempre = (temp-vmin)/(vmm+(vmm==0));
    %val_temp = vals(v);
    %vectre(:,:,v)=tempre/sqrt(val_temp+(val_temp==0));
end

end