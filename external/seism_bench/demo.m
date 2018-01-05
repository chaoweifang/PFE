% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 
%  Code obtained from:
%  https://imatge.upc.edu/web/resources/supervised-evaluation-image-segmentation
% ------------------------------------------------------------------------
% This file is part of the SEISM package presented in:
%    Jordi Pont-Tuset, Ferran Marques,
%    "Measures and Meta-Measures for the Supervised Evaluation of Image Segmentation,"
%    Computer Vision and Pattern Recognition (CVPR), 2013.
% If you use this code, please consider citing the paper.
% ------------------------------------------------------------------------ 
%
% Simple demo of how to use the package.
% For a realistic use, see "pr_curves.m".
%
% ------------------------------------------------------------------------ 

% List of measures to compute
measures = {%
            'fb'  ,... % Precision-recall for boundaries
            'fop' ,... % Precision-recall for objects and parts
            'fr'  ,... % Precision-recall for regions
            'voi' ,... % Variation of information
            'nvoi',... % Normalized variation of information
            'pri' ,... % Probabilistic Rand index
            'sc'  ,'ssc' ,... % Segmentation covering (two directions)
            'dhd' ,'sdhd',... % Directional Hamming distance (two directions)
            'bgm' ,... % Bipartite graph matching
            'vd'  ,... % Van Dongen
            'bce' ,... % Bidirectional consistency error
            'gce' ,... % Global consistency error
            'lce' ,... % Local consistency error
            };
        
% Compare a ground truth partition against the others
load(fullfile(root_dir, 'bsds500','ground_truth','120003.mat'))
partition = gt_seg{1};
gt        = gt_seg(2:3);

% Compute all available measures
home
disp('Computing measures...')
for ii=1:length(measures)
    tic
    result = eval_segm(partition, gt, measures{ii});
    disp([measures{ii} repmat(' ',1,6-length(measures{ii})) ': ' num2str(result) '   (' sprintf('%.3f',toc) 's.)'])
end
disp('Done!')

addpath('../PFE_MCG/full_new/external/bench/benchmarks/')
[ri,voi] = match_segmentations2(partition, gt);

matches = segmentationst(partition,gt);


tic
results = mex_eval_segm_pvc(partition,gt);
toc

regionsGT = [];
total_gt = 0;
gt1=cell(numel(gt),1);
nsegs = numel(gt);
for s = 1 : nsegs
    gt1{s} = double(gt{s});
    regionsTmp = regionprops(gt{s}, 'Area');
    regionsGT = [regionsGT; regionsTmp];
    total_gt = total_gt + max(gt{s}(:));
end


matchesSeg = max(matches, [], 2);
matchesGT = max(matches, [], 1);
cntP=0;
sumP=0;
cntR=0;
sumR=0;
regionsSeg = regionprops(partition, 'Area');
for r = 1 : numel(regionsSeg)
    cntP = cntP + regionsSeg(r).Area*matchesSeg(r);
    sumP = sumP + regionsSeg(r).Area;
end
cmgt = zeros(size(regionsGT));
for r = 1 : numel(regionsGT),
    cntR = cntR +  regionsGT(r).Area*matchesGT(r);
    cmgt(r) = regionsGT(r).Area*matchesGT(r);
    sumR = sumR + regionsGT(r).Area;
end
