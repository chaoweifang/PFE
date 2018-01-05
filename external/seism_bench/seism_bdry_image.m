function seism_bdry_image(inFile, gtFile, evFile, nthresh)
% Calculate region benchmarks for an image. Probabilistic Rand Index, Variation of
% Information and Segmentation Covering. 
%
% INPUT
%	inFile  : Can be one of the following:
%             - a collection of segmentations in a cell 'segs' stored in a mat file
%             - an ultrametric contour map in 'doubleSize' format, 'ucm2'
%               stored in a mat file with values in [0 1].
%
%	gtFile	:   File containing a cell of ground truth segmentations
%   evFile2, evFile3, evFile4  : Temporary output for this image.
%	nthresh:    Number of scales evaluated. If input is a cell of
%               'segs', nthresh is changed to the actual number of segmentations

if nargin<4, nthresh = 99; end

load(inFile); 
if exist('ucm2', 'var'),
    thresh = (1:nthresh)'/(nthresh+1);
elseif exist('segs', 'var')
    nthresh = numel(segs);
    thresh = (1:nthresh)';
else
    error('unexpected input in inFile');
end
%thresh=0.1;
load(gtFile);
if exist('groundTruth','var')
    nsegs = numel(groundTruth);
    if nsegs == 0,
        error(' bad gtFile !');
    end
    gt={};
    for i=1:nsegs
        gt = cat(1,gt,groundTruth{i}.Segmentation);
    end
elseif exist('LabelMap','var')
    gt{1} = LabelMap;
elseif exist('GTinst','var')
    gt{1} = GTinst.Segmentation+1;
else
    error('wrong input gt data or directory!');
end

measures = [];
for t=1:nthresh
    if exist('segs', 'var')
        seg = relabel(segs{t});
    else
        %labels2 = bwlabel(ucm <= thresh(t));
        %seg = labels2(2:2:end, 2:2:end);
        tmp = (bwlabel(ucm2'<=thresh(t),4))';
        seg = tmp(2:2:end,2:2:end);
    end
    [pre, rec] = fb(seg, gt);
    tpr = pre+rec;
    fpr = 2*pre*rec/(tpr+(tpr==0));
    measures = cat(1,measures,[fpr pre rec]);
end
%fop(value, fprec, frec) 
dlmwrite(evFile, measures);
end