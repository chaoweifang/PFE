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
using namespace std;
const int NM = 13;
const string STR_MEASURES[NM]={"fop","pri","fr","voi","nvoi","bce","gce","lce","sc","bgm","dhd","ssc","sdhd"};
enum CMeasure{
    M_FOP,
    M_PRI,
    M_FR,
    M_VOI,
    M_NVOI,
    M_BCE,
    M_GCE,
    M_LCE,
    M_SC,
    M_BGM,
    M_DHD,
    M_SSC,
    M_SDHD
};
void
mexFunction( int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[] )
{
    // Check number of inputs
    if (nrhs<3)
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
    
    bool measure_ids[NM];
    for(int i=0; i<NM; i++) 
    {
        measure_ids[i]=false;
    }

    for(int i=2; i<nrhs; i++)
    {
        std::string measure_id = "";
        if ( !mxIsChar(prhs[i]) )
            mexErrMsgTxt("Measure id must be a string.");
        else
            measure_id = mxArrayToString(prhs[i]);
        for(int j=0; j<NM; j++)
        {
            if(!strcmp(measure_id.c_str(),STR_MEASURES[j].c_str()))
            {
                measure_ids[j]=true;
            }
        }

    }
    //bool fop=false, pri=false, fr=false, voi=false,nvoi=false,bce=false,gce=false,lce=false,sc=false,bgm=false,dhd=false,ssc=false,sdhd=false;

    // Get measure name
    
    
    // Compute intersection matrices
    std::vector<MultiArray<uint64,2> > int_mats(n_gts);
    for(std::size_t ii=0; ii<n_gts; ++ii)
         int_mats[ii] = intersection_matrix(part, gts[ii]);;
    
    // Vector to store output
    std::vector<double> results;

    /***********************************************/
    /*        Objects and Parts F-measure          */
    /***********************************************/
    if (measure_ids[M_FOP])
    {
        float64 gamma_obj  = 0.9;
        float64 gamma_part = 0.25;
        float64 beta       = 0.1;
        ObjectsAndPartsF opf(gamma_obj, gamma_part, beta);
        opf.calculate(int_mats);
        results.push_back(opf.f_measure());
        results.push_back(opf.precision());
        results.push_back(opf.recall());
    }
    
    
    /*****************************************************/
    /* Probabilistic Rand Index (PRI): Mean of RI to GTs */
    /*****************************************************/
    if (measure_ids[M_PRI])
    {
        double pri = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            PairsOfPixels pairs;
            pairs.calculate(int_mats[ii]);
            pri += (double)pairs.rand_index();
        }
        pri /= (double)int_mats.size();
        results.push_back(pri);
    }
    
    /***********************************************/
    /*        Precision-recall for regions         */
    /***********************************************/
    if (measure_ids[M_FR])
    {
        float64 prr = 0; float64 rer = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            PairsOfPixels pairs;
            pairs.calculate(int_mats[ii]);
            prr += pairs.precision();
            rer += pairs.recall();
        }
        prr /= (float64)n_gts;
        rer /= (float64)n_gts;
        float64 fr = 0;
        if (prr+rer>0)
            fr = 2*prr*rer/(prr+rer);
        else
            fr = 0;
        results.push_back(fr); 
        results.push_back(prr);
        results.push_back(rer);
    }    
    
    /***********************************************/
    /*          Variation of Information           */
    /***********************************************/
    if (measure_ids[M_VOI])
    {
        float64 voi = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            voi += (float64)variation_of_information(int_mats[ii]);
        }
        voi /= (float64)n_gts;
        results.push_back(voi); 
    }
    
    /***********************************************/
    /*    Normalized Variation of Information      */
    /***********************************************/
    if (measure_ids[M_NVOI])
    {
        float64 nvoi = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            uint32 num_regs_1 = int_mats[ii].dims(0);
            uint32 num_regs_2 = int_mats[ii].dims(1);

            float64 tmp = (float64)variation_of_information(int_mats[ii]);

            // Normalized VoI
            if (std::max(num_regs_1,num_regs_2)==1)
                nvoi += 1;
            else
                nvoi += 1-(tmp/(2*log(std::max(num_regs_1,num_regs_2))/log(2)));
        }
        nvoi /= (float64)n_gts;
        results.push_back(nvoi);
    }
    
    /***********************************************/
    /*       Bidirectional Consistency Error       */
    /***********************************************/
    if (measure_ids[M_BCE])
    {
        float64 bce = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            bce += (float64)bidirectional_consistency_error(int_mats[ii]);
        }
        bce /= (float64)n_gts;
        bce = 1-bce;  // Make it a similarity measure (1->Best)
        results.push_back(bce); 
    }
    
    /***********************************************/
    /*         Global Consistency Error            */
    /***********************************************/
    if (measure_ids[M_GCE])
    {
        float64 gce = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            gce += (float64)global_consistency_error(int_mats[ii]);
        }
        gce /= (float64)n_gts;
        gce = 1-gce;  // Make it a similarity measure (1->Best)
        results.push_back(gce); 
    }    
    
    /***********************************************/
    /*         Local Consistency Error             */
    /***********************************************/
    if (measure_ids[M_LCE])
    {
        float64 lce = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            lce += (float64)local_consistency_error(int_mats[ii]);
        }
        lce /= (float64)n_gts;
        lce = 1-lce;  // Make it a similarity measure (1->Best)
        results.push_back(lce); 
    }   
    
    /***********************************************/
    /*          Segmentation Covering              */
    /***********************************************/
    if (measure_ids[M_SC])
    {
        float64 sc = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            sc += (float64)segmentation_covering(int_mats[ii]);
        }
        sc /= (float64)n_gts;
        results.push_back(sc); 
    }   
    
    /***********************************************/
    /*       Bipartite Graph Matching (BGM)        */
    /***********************************************/
    if (measure_ids[M_BGM])
    {
        float64 bgm = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            bgm += (float64)bipartite_graph_matching(part, gts[ii]);
        }
        bgm /= (float64)n_gts;  // Mean over GTs
        bgm /= (float64)(part.dims()[0]*part.dims()[1]);  // Normalize
        bgm = 1-bgm;  // Make it a similarity measure (1->Best)
        results.push_back(bgm);
    }
    
    /***********************************************/
    /*    Directional Hamming Distance (DHD)       */
    /***********************************************/
    if (measure_ids[M_DHD])
    {
        float64 dhd = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            dhd += (float64)directional_hamming_distance(int_mats[ii]);
        }
        dhd /= (float64)n_gts;
        dhd /= (float64)(part.dims()[0]*part.dims()[1]);  // Normalize
        dhd = 1-dhd;  // Make it a similarity measure (1->Best)
        results.push_back(dhd); 
    }
    
    /***********************************************/
    /*     Swapped Segmentation Covering (SSC)     */
    /***********************************************/
    if (measure_ids[M_SSC])
    {
        float64 sc = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            /* Transpose int_mat */
            MultiArray<uint64,2> transp_int(int_mats[ii].dims()[1],int_mats[ii].dims()[0]);
            for(std::size_t xx=0; xx<int_mats[ii].dims()[0]; xx++)
                for(std::size_t yy=0; yy<int_mats[ii].dims()[1]; yy++)
                    transp_int[yy][xx] = int_mats[ii][xx][yy];
            
            sc += (float64)segmentation_covering(transp_int);
        }
        sc /= (float64)n_gts;
        results.push_back(sc);
    }
    
    /***********************************************/
    /* Swapped Directional Hamming Distance (SDHD) */
    /***********************************************/
    if (measure_ids[M_SDHD])
    {
        float64 dhd = 0;
        for(std::size_t ii=0; ii<n_gts; ++ii)
        {
            /* Transpose int_mat */
            MultiArray<uint64,2> transp_int(int_mats[ii].dims()[1],int_mats[ii].dims()[0]);
            for(std::size_t xx=0; xx<int_mats[ii].dims()[0]; xx++)
                for(std::size_t yy=0; yy<int_mats[ii].dims()[1]; yy++)
                    transp_int[yy][xx] = int_mats[ii][xx][yy];
            
            dhd += (float64)directional_hamming_distance(transp_int);
        }
        dhd /= (float64)n_gts;
        dhd /= (float64)(part.dims()[0]*part.dims()[1]);  // Normalize
        dhd = 1-dhd;  // Make it a similarity measure (1->Best)
        results.push_back(dhd); 
    }
    
    // Fill output
    plhs[0]  =  mxCreateDoubleMatrix(1, results.size(), mxREAL);
    MatlabMultiArray<double> resul(plhs[0]);
    for(std::size_t ii=0; ii<results.size(); ++ii)
    {
        resul[0][ii] = results[ii];
    }
}
