function im_vec=showEmbedding(vect,w,h)
pad=10;
chs=size(vect,2);
im_vec=ones(h,w*chs+pad*(chs-1));
for i=1:chs
    vec=vect(:,i);
    vec=(vec-min(vec(:)))/(max(vec(:))-min(vec(:)));
    vec=reshape(vec,h,w);
    im_vec(1:h,(i-1)*(w+pad)+1:(i-1)*(w+pad)+w)=vec;
end
end

