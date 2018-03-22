/*
*Piecewise flat embedding for image segmentation
*GMM Estimation Using multithreads
*@Chaowei FANG
*12/3/2015
*/

// includes, system
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h> // for clock(), clock_t, CLOCKS_PER_SEC

#include "gaussian.h"
#include "invert_matrix.h"
#include "gaussian_kernel.h"



////////////////////////////////////////////////////////////////////////////////
// Program mexFunction
////////////////////////////////////////////////////////////////////////////////
void mexFunction(
	int nargout,
	mxArray *out[],
	int nargin,
	const mxArray *in[]
	)
{
	int num_clusters;   //number of cluster centers
	int num_dimensions; //dimension of input data
	int num_events;     //number of input data
	int num_threads;    //number of threads to compute mstep
	//io_timer = clock();
	if (nargin != 3 && nargin!=5)
		mexErrMsgTxt("3 or 5 parameters are required: input data and number of clusters!");
	if(!mxIsDouble(in[0]) || !mxIsDouble(in[1]) || !mxIsDouble(in[2]))
		mexErrMsgTxt("input data and number of clusters should be double type!");
	num_dimensions = (int)mxGetM(in[0]);
	num_events = (int)mxGetN(in[0]);
	num_clusters = (int)mxGetScalar(in[1]);
	num_threads = (int)mxGetScalar(in[2]);
	int mstep_threads = num_threads<num_clusters? num_threads:num_clusters;
	int estep1_threads = mstep_threads;
	int estep2_threads = num_threads<num_events? num_threads:num_events;

	double *inicenters, *inicov;  //initialized centers and covariance
	int num_c;           
	const size_t* dims_covar;
	int ndims;
	if(nargin > 3)
	{
		if((int)mxGetM(in[3]) != num_dimensions || (int)mxGetN(in[3]) != num_clusters || !mxIsDouble(in[3]))
			mexErrMsgTxt("Inilialized centers are required to be num_dimensions*num_clusters double matrix");
		inicenters = mxGetPr(in[3]);
	}
	if(nargin==5)
	{
		ndims = (int)mxGetNumberOfDimensions(in[4]);
		if(ndims != 3 || !mxIsDouble(in[4]))
			mexErrMsgTxt("Initialized covariance matrix is required to be double with 3 channels");
		dims_covar = mxGetDimensions(in[4]);
		if(dims_covar[0]!=dims_covar[1] || (int)dims_covar[0]!=num_dimensions || (int)dims_covar[2]!=num_clusters)
			mexErrMsgTxt("Wrong input covariance matrix");
		inicov = mxGetPr(in[4]);
	}


	// Transpose the event data (allows coalesced access pattern in E-step kernel)
	// This has consecutive values being from the same dimension of the data 
	// (num_dimensions by num_events matrix)
	double* fcs_data_by_event = mxGetPr(in[0]);
	double* fcs_data_by_dimension = (double*)malloc(sizeof(double)*num_events*num_dimensions);

	for (int e = 0; e<num_events; e++) {
		for (int d = 0; d<num_dimensions; d++) {
			fcs_data_by_dimension[d*num_events + e] =(double) fcs_data_by_event[e*num_dimensions + d];
		}
	}

	out[0] = mxCreateNumericMatrix(num_dimensions,num_clusters,mxDOUBLE_CLASS,mxREAL); //means
	mwSize cov_dims[3]; cov_dims[0]=num_dimensions; cov_dims[1]=num_dimensions; cov_dims[2]=num_clusters;               //
	out[1] = mxCreateNumericArray(3,cov_dims,mxDOUBLE_CLASS,mxREAL);             //covariances matrixes
	out[2] = mxCreateNumericMatrix(1,num_clusters,mxDOUBLE_CLASS,mxREAL);            //weight vector
	out[3] = mxCreateNumericMatrix(num_events,num_clusters,mxDOUBLE_CLASS,mxREAL);   //
	out[4] = mxCreateDoubleScalar(0);
	double* ptr_lh=mxGetPr(out[4]);

	// Setup the cluster data structures on host
	clusters_t clusters;
	clusters.N = (double*)malloc(sizeof(double)*num_clusters);
	clusters.pi = mxGetPr(out[2]);
	clusters.constant = (double*)malloc(sizeof(double)*num_clusters);
	clusters.avgvar = (double*)malloc(sizeof(double)*num_clusters);
	clusters.means = mxGetPr(out[0]);
	clusters.R = mxGetPr(out[1]);
	clusters.Rinv = (double*)malloc(sizeof(double)*num_dimensions*num_dimensions*num_clusters);
	clusters.memberships = mxGetPr(out[3]);
	if (!clusters.means || !clusters.R || !clusters.Rinv || !clusters.memberships) {
		mexErrMsgTxt("ERROR: Could not allocate memory for clusters.\n");
	}

	double rissanen;

	//////////////// Initialization done, starting kernels //////////////// 
	//DEBUG("Invoking seed_clusters kernel.\n");
	//fflush(stdout);

	// seed_clusters sets initial pi values, 
	// finds the means / covariances and copies it to all the clusters
	// TODO: Does it make any sense to use multiple blocks for this?
	//seed_start = clock();
	if(nargin != 5)
	{
		seed_clusters(fcs_data_by_event, &clusters, num_dimensions, num_clusters, num_events);
	}
	else
	{
		memcpy(clusters.means, inicenters, num_clusters*num_dimensions*sizeof(double));
		memcpy(clusters.R,inicov,num_clusters*num_dimensions*num_dimensions*sizeof(double));
		for(int m=0; m < num_clusters; m++) {
			clusters.N[m] = (double) num_events / (double) num_clusters;
        	clusters.pi[m] = 1.0 / num_clusters;
		}
	}

	//DEBUG("Invoking constants kernel.\n");
	// Computes the R matrix inverses, and the gaussian constant
	//constants_kernel<<<num_clusters, num_threads>>>(d_clusters,num_clusters,num_dimensions);
	constants(&clusters, num_clusters, num_dimensions);

	// Calculate an epsilon value
	//int ndata_points = num_events*num_dimensions;
	double epsilon = (1 + num_dimensions + 0.5f*(num_dimensions + 1)*num_dimensions)*(double)log((double)num_events*num_dimensions)*0.01f;
	double likelihood, old_likelihood;
	int iters;
	epsilon = 1e-6f;

	/*************** EM ALGORITHM *****************************/

	// do initial regrouping
	// Regrouping means calculate a cluster membership probability
	// for each event and each cluster. Each event is independent,
	// so the events are distributed to different blocks 
	// (and hence different multiprocessors)
	mt_estep1(fcs_data_by_dimension, &clusters, num_dimensions, num_clusters, num_events, estep1_threads);
	mt_estep2(fcs_data_by_dimension, &clusters, &likelihood, num_dimensions, num_clusters, num_events, estep2_threads);
	
	double change = epsilon * 2;
	iters = 0;     
	while (iters < MIN_ITERS || (fabs(change) > epsilon && iters < MAX_ITERS)) {
		old_likelihood = likelihood;
		mt_mstep(fcs_data_by_dimension, &clusters, num_dimensions, num_clusters, num_events,mstep_threads);
		constants(&clusters, num_clusters, num_dimensions);
		mt_estep1(fcs_data_by_dimension, &clusters, num_dimensions, num_clusters, num_events, estep1_threads);
		mt_estep2(fcs_data_by_dimension, &clusters, &likelihood, num_dimensions, num_clusters, num_events, estep2_threads);
		change = likelihood - old_likelihood;


		iters++;

	}
	
	// Calculate Rissanen Score
	rissanen = -likelihood + 0.5*(num_clusters*(1 + num_dimensions + 0.5*(num_dimensions + 1)*num_dimensions) - 1)*log(num_events*num_dimensions);
	*ptr_lh = rissanen;
	
	free(clusters.N);
	free(clusters.constant);
	free(clusters.avgvar);
	free(clusters.Rinv);
	free(fcs_data_by_dimension);
}

