function [seg_best,fop_best]=seism_reg_simg(data, gtFile, nthresh,itype)
% Calculate region benchmarks for an image. Probabilistic Rand Index, Variation of
% Information and Segmentation Covering. 
%
% INPUT
%	data  : Can be one of the following:
%             - a collection of segmentations in a cell 'segs' stored in a mat file
%             - an ultrametric contour map in 'doubleSize' format, 'ucm2'
%               stored in a mat file with values in [0 1].
%
%	gtFile	:   File containing a cell of ground truth segmentations
%	nthresh:    Number of scales evaluated. If input is a cell of
%               'segs', nthresh is changed to the actual number of segmentations
%   itype:   indicate type of data, 'ucm2'/'segs'
if ~exist('nthresh','var')||isempty(nthresh), nthresh = 100; end
if strcmp(itype,'ucm2')
    thresh = (1:nthresh)'/(nthresh+1); 
else if strcmp(itype,'segs')
        nthresh = numel(data);
        thresh = (1:nthresh)';
    else
        error('unexpected input type!!!')
    end
end

load(gtFile);
nsegs = numel(groundTruth);
if nsegs == 0,
    error(' bad gtFile !');
end
gt={};
for i=1:nsegs
    gt = cat(1,gt,groundTruth{i}.Segmentation);
end
%measures = [];
seg_best = [];
fop_best = 0;
for t=1:nthresh
    if exist('segs', 'var')
        seg = data{t};
    else
        %labels2 = bwlabel(ucm <= thresh(t));
        %seg = labels2(2:2:end, 2:2:end);
        tmp = (bwlabel(data'<=thresh(t),4))';
        seg = tmp(2:2:end,2:2:end);
    end
    measure_tmp = eval_segm(seg, gt,'fop');
    if measure_tmp(1) > fop_best
        fop_best = measure_tmp(1);
        seg_best = seg;
    end
    %measures = cat(1,measures,measure_tmp);
end
%fop(value, fprec, frec) 
%dlmwrite(evFile, measures);
end