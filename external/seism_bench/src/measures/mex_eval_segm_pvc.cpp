// ------------------------------------------------------------------------ 
//  Copyright (C)
//  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
// 
//  Jordi Pont-Tuset <jordi.pont@upc.edu>
//  June 2013
// ------------------------------------------------------------------------ 
//  Code available at:
//  https://imatge.upc.edu/web/resources/supervised-evaluation-image-segmentation
// ------------------------------------------------------------------------
//  This file is part of the SEISM package presented in:
//    Jordi Pont-Tuset, Ferran Marques,
//    "Measures and Meta-Measures for the Supervised Evaluation of Image Segmentation,"
//    Computer Vision and Pattern Recognition (CVPR), 2013.
//  If you use this code, please consider citing the paper.
// ------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <aux/intersection_matrix.hpp>
#include "mex.h"

#include <measures/fop.hpp>
#include <measures/pri_fr.hpp>
#include <measures/bgm.hpp>
#include <measures/voi.hpp>
#include <measures/bce_gce_lce.hpp>
#include <measures/sc.hpp>
#include <measures/dhd.hpp>

void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{
    // Check number of inputs
    if (nrhs<2)
        mexErrMsgTxt("At least two arguments needed");

    // Check type of partition
    if (mxGetClassID(prhs[0])!=mxUINT32_CLASS)
        mexErrMsgTxt("The partition should be in uint32\n");
    
    // Partition as multi_array
    MultiArray<uint32,2> part;
    part = ConstMatlabMultiArray<uint32>(prhs[0]);
       
    // Check type of ground truth
    if (mxGetClassID(prhs[1])!=mxCELL_CLASS)
        mexErrMsgTxt("Ground truth should be a cell\n");
    std::size_t n_gts = mxGetNumberOfElements(prhs[1]);

    // Ground truths as multi_array
    std::vector<MultiArray<uint32,2> > gts(n_gts);
    for(std::size_t ii=0; ii<n_gts; ++ii)
    {
        if (mxGetClassID(mxGetCell(prhs[1],ii))!=mxUINT32_CLASS)
            mexErrMsgTxt("The ground truths should be in uint32\n");
        
        ConstMatlabMultiArray<uint32> tmp(mxGetCell(prhs[1],ii));
        gts[ii] = tmp;
    }
    
    // Get measure name
    //std::string measure_id = "fop";
    //if ( mxIsChar(prhs[2]) != 1)
     //   mexErrMsgTxt("Measure id must be a string.");
    //else
     //   measure_id = mxArrayToString(prhs[2]);
    
    // Compute intersection matrices
    std::vector<MultiArray<uint64,2> > int_mats(n_gts);
    for(std::size_t ii=0; ii<n_gts; ++ii)
         int_mats[ii] = intersection_matrix(part, gts[ii]);;
    
    // Vector to store output
    std::vector<double> results;

    
    
    /*****************************************************/
    /* Probabilistic Rand Index (PRI): Mean of RI to GTs */
    /*****************************************************/
        double pri = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            PairsOfPixels pairs;
            pairs.calculate(int_mats[ii]);
            pri += (double)pairs.rand_index();
        }
        pri /= (double)int_mats.size();
        results.push_back(pri);
    
    /***********************************************/
    /*          Variation of Information           */
    /***********************************************/
        float64 voi = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            voi += (float64)variation_of_information(int_mats[ii]);
        }
        voi /= (float64)n_gts;
        results.push_back(voi); 
    
    /***********************************************/
    /*          Segmentation Covering              */
    /***********************************************/
       /*uint32 num_reg_1 = int_mats[0].dims()[0];
        for(size_t i=0; i<n_gts; i++)
    {
        const MultiArray<uint64,2> intersect_matrixt = int_mats[i];
        uint32 num_reg_2t = int_mats[i].dims()[1];
        printf("%d \n",num_reg_2t);
        for(std::size_t jj = 0; jj < num_reg_2t; jj++)
        {
            for(std::size_t ii = 0; ii < num_reg_1; ii++)
            {
                 printf("i=%d jj=%d ii=%d inter=%d \n",i,jj,ii,(int)intersect_matrixt[ii][jj]);
            }
        }
    }*/
        //const MultiArray<uint64,2> intersect_matrixt = pintersect_matrix[i];
        std::vector<float64> result_sc;
        segmentation_covering_count(int_mats,result_sc);
        results.insert(results.end(),result_sc.begin(),result_sc.end());
       /* float64 sc = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            sc += (float64)segmentation_covering(int_mats[ii]);
        }
        sc /= (float64)n_gts;
        results.push_back(sc); 

    /***********************************************/
    /*     Swapped Segmentation Covering (SSC)     */
    /***********************************************/
    /*    float64 sc = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            /* Transpose int_mat */
            /*MultiArray<uint64,2> transp_int(int_mats[ii].dims()[1],int_mats[ii].dims()[0]);
            for(std::size_t xx=0; xx<int_mats[ii].dims()[0]; xx++)
                for(std::size_t yy=0; yy<int_mats[ii].dims()[1]; yy++)
                    transp_int[yy][xx] = int_mats[ii][xx][yy];
            
            sc += (float64)segmentation_covering(transp_int);
        }
        sc /= (float64)n_gts;
        results.push_back(sc);*/
    
    // Fill output
    plhs[0]  =  mxCreateDoubleMatrix(1, results.size(), mxREAL);
    MatlabMultiArray<double> resul(plhs[0]);
    for(std::size_t ii=0; ii<results.size(); ++ii)
    {
        resul[0][ii] = results[ii];
    }
}
