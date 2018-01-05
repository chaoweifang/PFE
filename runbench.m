addpath('./datasets')
addpath(genpath('./external/seism_bench/'))
database='VOCContext';
dset='test';
method = 'swPFE0.8';
resroot='/home2/cwfang/Results';

indir = fullfile(resroot,database,method);
outdir =  [indir '/seism_test'];
seism_region(database,dset,indir,outdir)
seism_plot_reg(outdir)
seism_bdry(database,dset,indir,outdir);
seism_plot_reg(outdir)
