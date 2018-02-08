# Piecewise Flat Embedding for Image Segmentation


![](https://ws4.sinaimg.cn/large/006tNc79gy1fn3p0c3wd6j31hw0iadro.jpg)

Piecewise Flat Embedding (PFE) is a new nonlinear embedding for image segmentation. Based on the theory of sparse signal recovery, piecewise flat embedding attempts to recover a piecewise constant image representation with sparse region boundaries and sparse cluster value scattering. The resultant piecewise flat embedding exhibits interesting properties such as suppressing slowly varying signals, and offers an image representation with higher region identifiability which is desirable for image segmentation or high-level semantic analysis tasks.

Some examples of embedding results:
![](https://ws3.sinaimg.cn/large/006tNc79gy1fn3p73b8cfj30s50f9tl2.jpg)

### Prerequisites

A boundary detector is required to generate boundary probability map which is used to build the affinity matrix in PFE. There are piles of choices for boundary detection such as [HED](https://github.com/s9xie/hed), [DB](http://www0.cs.ucl.ac.uk/staff/I.Kokkinos/deepboundaries.html) and [COB](http://www.vision.ee.ethz.ch/~cvlsegmentation/cob/).

### Usage

1. run_im2ucm_all.m: generate and save segmentation results for all images the chosen dataset;
2. run_bench.m: evaluate segmentation results using boundary (Fb) and segmenation (Fop, Covering, PRI and VI) metrics;
3. different versions of PFE can be set in get_parameters.m. Set p=1 to use the L1-regularized PFE in demo.m.


### Precomputed Results

Contour driven segmentations of methods integrating swPFE0.8 precomputed on [Standford Background Dataset(SBD)](http://dags.stanford.edu/projects/scenedataset.html), [MSRC](http://riemenschneider.hayko.at/vision/dataset/task.php?did=36), and testing subset of [BSDS500](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html) and [PASCALContext](http://www.cs.stanford.edu/~roozbeh/pascal-context/), can be downloaded in the following links. 

1. Use COB as local bounadry: [BSDS500](https://drive.google.com/open?id=1gdMyp07kcR9OpbSSkEkV7MtpdP2qXrSL), [SBD](https://drive.google.com/open?id=1X1uJEfpeyboxOIwtXRsY-MGb6kDabTDw), [MSRC](https://drive.google.com/open?id=1C0j_658lJ_MGuqadPBM2AyRJKaD0iA0m), [PASCALContext](https://drive.google.com/open?id=19GasnZOi_hLW8INbV4WNqlZHtgvCVflp);
2. Use DB as local boundary: [BSDS500](https://drive.google.com/open?id=16pwlZOmGFSTxzLqLneeDWK_Xowf4QTe_), [SBD](https://drive.google.com/open?id=1_acRQiucHu4qgstXnySoAO4h6IvutWKN), [MSRC](https://drive.google.com/open?id=1YoLhmOZ2rLZN9C2Sy2cW8Wbfx1IiH_VA), [PASCALContext](https://drive.google.com/open?id=1cGG_BMqXeAJY1B6b5UlXR9fpgrl2Lvp_).

### Citation

If you use this code, please citing the following papers:

@inproceedings{yu2015piecewise,title={Piecewise flat embedding for image segmentation},
  author={Yu, Yizhou and Fang, Chaowei and Liao, Zicheng},
  booktitle={Proceedings of the IEEE International Conference on Computer Vision},
  pages={1368--1376},
  year={2015}
}

Enjoy!
