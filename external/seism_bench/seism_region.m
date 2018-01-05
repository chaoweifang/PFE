function seism_region(database,dset,inDir,outDir,nthresh)
% Run region benchmarks on dataset: Probabilistic Rand Index, Variation of
% Information and Segmentation Covering. 
%
% INPUT
%   imgDir: folder containing original images
%   gtDir:  folder containing ground truth data.
%   inDir:  folder containing segmentation results for all the images in imgDir. 
%           Format can be one of the following:
%             - a collection of segmentations in a cell 'segs' stored in a mat file
%             - an ultrametric contour map in 'doubleSize' format, 'ucm2' stored in a mat file with values in [0 1].
%   outDir: folder where evaluation results will be stored
%	nthresh	: Number of points in precision/recall curve.

%imgDir ='~/MATLAB/BSDS500/data/images/test';%'../../Val/images'; 
%gtDir =  '~/MATLAB/BSDS500/data/groundTruth/test';
%inDir = '~/MATLAB/Data/MCG-BSDS500-all-ucm';
%outDir = fullfile(inDir,'test_seism');
if ~exist(outDir,'dir'), mkdir(outDir); end
if exist(fullfile(outDir,'eval_fop_pr.txt'),'file')
    disp('already evaluated!');
    return
end
if nargin<5, nthresh=99; end

root= database_root_dir(database);
gtDir = fullfile(root,'gt');
iids = database_ids(database, dset);
%hpool = parpool(16);
parfor i = 1 : numel(iids)
    str_imgid = iids{i};
    evFile = fullfile(outDir, strcat(str_imgid, '_ev.txt'));
    if ~isempty(dir(evFile)), continue; end
    inFile = fullfile(inDir, strcat(str_imgid, '.mat'));
    gtFile = fullfile(gtDir, strcat(str_imgid, '.mat'));
    if ~exist(inFile,'file') || ~exist(gtFile,'file')
        error(['input or gt file not exists for ' str_imgid])
    end
    seism_reg_image(inFile, gtFile, evFile, nthresh);
end

%get thresh array
thresh = [];
if numel(iids)>0
    load(fullfile(inDir, strcat(iids{1}, '.mat')));
    if exist('ucm2', 'var'),
        thresh = (1:nthresh)'/(nthresh+1);
        clear ucm2
    elseif exist('segs', 'var')
        nthresh = numel(segs);
        thresh = (1:nthresh)';
        clear segs
    else
        error('unexpected input in inFile');
    end
end

%collect benchmarks
collect_seism_reg(outDir,thresh)
delete(fullfile(outDir,'*_ev.txt'));
end