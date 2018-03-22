/*
 * Parameters file for gaussian mixture model based clustering application
 */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

// Maxinum number of threads per block is 512, so that limits us to 512 clusters
// Probably will run out of memory and make the computation intractable far before 512 clusters though
#define MAX_CLUSTERS 512
#define PI  3.1415926535897931
#define	NUM_BLOCKS 30
#define NUM_THREADS 512 // Must be power of 2 due to butterfly sum reductions
#define NUM_DIMENSIONS 32

#define COVARIANCE_DYNAMIC_RANGE 1E6

// if 0, uses random, else picks events uniformly distributed in data set
#define UNIFORM_SEED 1

// Which GPU to use, if more than 1
#define DEVICE 1

// Using only diagonal covariance matrix, thus all dimensions are considered independent
#define DIAG_ONLY 0

// Maximum number of iterations for the EM convergence loop
#define MIN_ITERS 0
#define MAX_ITERS 100

typedef struct 
{
    // Key for array lengths
    //  N = number of events
    //  M = number of clusters
    //  D = number of dimensions
    double* N;        // expected # of pixels in cluster: [M]
    double* pi;       // probability of cluster in GMM: [M]
    double* constant; // Normalizing constant [M]
    double* avgvar;    // average variance [M]
    double* means;   // Spectral mean for the cluster: [M*D]
    double* R;      // Covariance matrix: [M*D*D]
    double* Rinv;   // Inverse of covariance matrix: [M*D*D]
    double* memberships; // Fuzzy memberships: [N*M]
} clusters_t;

typedef struct 
{
    double* data;               // input data
    double* likelihood;
    clusters_t* clusters;       // cluster data
    int D;                      // number of dimensions
    int M;                      // number of clusters
    int N;                      // number of events
    int td_id;                  // thread id
    int len;
} thread_t;

#endif

