function [ucm2,sPb,vecs,Yor,t]= img2ucm(I,vecs,Yor,pb,param_s,param_multi,param_ini)
%initialization
tic
if isempty(Yor)
    img_t=[];
    img = im2double(I);
    for i=1:size(img,3)
       img_t = cat(3,img_t,img(:,:,i)');
    end
    if strcmp(param_ini.flag,'g')
        Yor = IniMtGMM(img_t,param_ini);
    else
        if strcmp(param_ini.flag,'h')
            Yor = IniGBGaussian(img_t,param_ini);
        else
            if strcmp(param_ini.flag,'s')
                Yor = IniWSC(img_t,param_ini);
            else
                error('Error initialization type!');
            end
        end
    end
end
t_ini=toc;
tic
if isempty(vecs)
    O = edgeOrient(pb,param_s.or1);
    e_nms = edgesNmsMex(pb,O,1,5,1.01,4);
    lpb = e_nms * param_s.mult_Pb;
    vecs = pfe(lpb, Yor, param_s,param_ini.chs);
end
t_pfe=toc;
tic
[tx,ty]=size(pb);
if strcmp(param_s.wtype,'o') || (strcmp(param_s.wtype,'w') && param_s.p~=1)
    sPb = cvtVec2Edge(vecs.vect,tx,ty);
elseif strcmp(param_s.wtype,'w')
    sPb = cvtVec2EdgeRe(vecs,tx,ty);
end
sPb = sPb/param_s.sat_sPb;    
gPb = (double(pb)+sPb)/2;
On = edgeOrient(gPb,param_s.or2);
[owt2n, superpixels] = contours2OWT(gPb,On);

% ultrametric contour map with mean pb.
ucm = double(ucm_mean_pb( owt2n, superpixels));
ucm2 = apply_sigmoid(ucm,param_multi.thr,param_multi.fq); %ucm2multi(ucm,param_multi);
t_ucm=toc;
t_sum=t_ini+t_pfe+t_ucm;
t=[t_sum,t_ini,t_pfe,t_ucm];
end