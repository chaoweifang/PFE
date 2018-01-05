function ucm2=edg2ucm(pb,thr,fq)
O = edgeOrient(pb,3);
[owt2, superpixels] = contours2OWT(pb, O);
ucm = double(ucm_mean_pb( owt2, superpixels));
ucm2 = apply_sigmoid(ucm,thr,fq);