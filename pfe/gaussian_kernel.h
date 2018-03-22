#ifndef GAUSSIAN_KERNEL_H
#define GAUSSIAN_KERNEL_H
#include <algorithm>
#include <cfloat> 
#include "pthread.h"
//#include <ctgmath>
void seed_clusters(double *data, clusters_t* clusters, int D, int M, int N) {
    double* variances = (double*) malloc(sizeof(double)*D);
    double* means = (double*) malloc(sizeof(double)*D);

    // Compute means
    for(int d=0; d < D; d++) {
        means[d] = 0.0;
        for(int n=0; n < N; n++) {
            means[d] += data[n*D+d];
        }
        means[d] /= (double) N;
    }

    // Compute variance of each dimension
    for(int d=0; d < D; d++) {
        variances[d] = 0.0;
        for(int n=0; n < N; n++) {
            variances[d] += data[n*D+d]*data[n*D+d];
        }
        variances[d] /= (double) N;
        variances[d] -= means[d]*means[d];
    }

    // Average variance
    double avgvar = 0.0;
    for(int d=0; d < D; d++) {
        avgvar += variances[d];
    }
    avgvar /= (double) D;

    // Initialization for random seeding and uniform seeding    
    double fraction;
    //int seed;
    if(M > 1) {
        fraction = (N-1.0f)/(M-1.0f);
    } else {
        fraction = 0.0;
    }
    //srand(seed);
    srand(clock());

    for(int m=0; m < M; m++) {
        clusters->N[m] = (double) N / (double) M;
        clusters->pi[m] = 1.0f / (double) M;
        clusters->avgvar[m] = avgvar / COVARIANCE_DYNAMIC_RANGE;

        //DEBUG("N: %.2f\tPi: %.2f\tAvgvar: %e\n",clusters->N[m],clusters->pi[m],clusters->avgvar[m]);

        // Choose cluster centers
        //DEBUG("Means: ");
        #if UNIFORM_SEED
            for(int d=0; d < D; d++) {
                clusters->means[m*D+d] = data[((int)(m*fraction))*D+d];
                //DEBUG("%.2f ",clusters->means[m*D+d]);
            }
        #else
            seed = rand() % N;
            //DEBUG("Cluster %d seed = event #%d\n",m,seed);
            for(int d=0; d < D; d++) {
                clusters->means[m*D+d] = data[seed*D+d];
                //DEBUG("%.2f ",clusters->means[m*D+d]);
            }
        #endif
        //DEBUG("\n");

        // Set covariances to identity matrices
        for(int i=0; i < D; i++) {
            for(int j=0; j < D; j++) {
                if(i == j) {
                    clusters->R[m*D*D+i*D+j] = 1.0f;
                } else {
                    clusters->R[m*D*D+i*D+j] = 0.0f;
                }
            }
        }
        
        //DEBUG("R:\n");
        //for(int d=0; d < D; d++) {
        //    for(int e=0; e < D; e++) 
                //DEBUG("%.2f ",clusters->R[m*D*D+d*D+e]);
            //DEBUG("\n");
        //}
       // DEBUG("\n");

    }

    free(variances);
    free(means);
}

void constants(clusters_t* clusters, int M, int D) {
    double log_determinant;
    double* matrix = (double*) malloc(sizeof(double)*D*D);

    double sum = 0.0;
    for(int m=0; m < M; m++) {
        // Invert covariance matrix
        memcpy(matrix,&(clusters->R[m*D*D]),sizeof(double)*D*D);
        invert_cpu(matrix,D,&log_determinant);
        memcpy(&(clusters->Rinv[m*D*D]),matrix,sizeof(double)*D*D);

        // Compute constant
        clusters->constant[m] = -D*0.5*log(2.0*PI) - 0.5*log_determinant;
        //DEBUG("Cluster %d constant: %e\n",m,clusters->constant[m]);

        // Sum for calculating pi values
        sum += clusters->N[m];
    }

    // Compute pi values
    for(int m=0; m < M; m++) {
        clusters->pi[m] = clusters->N[m] / sum;
    }
    
    free(matrix);
}

void estep1(double* data, clusters_t* clusters, int D, int Mfirst, int Mlast, int N) {
    //clock_t start,finish;
    // Compute likelihood for every data point in each cluster
    double like;
    double* means;
    double* Rinv;
    int D_sq = D*D;
    int R_chn_idx = Mfirst*D_sq;
    int mem_col_idx = Mfirst*N;
    int mea_col_idx = Mfirst*D;
    //start = clock();
    for(int m=Mfirst; m < Mlast; m++,R_chn_idx+=D_sq,mem_col_idx+=N,mea_col_idx+=D) {
        means = (double*) &(clusters->means[mea_col_idx]);
        Rinv = (double*) &(clusters->Rinv[R_chn_idx]);
        
        for(int n=0; n < N; n++) {
            like = 0.0;
            #if DIAG_ONLY
                for(int i=0; i < D; i++) {
                    like += (data[i*N+n]-means[i])*(data[i*N+n]-means[i])*Rinv[i*D+i];
                }
            #else
                for(int i=0; i < D; i++) {
                    for(int j=0; j < D; j++) {
                        like += (data[i*N+n]-means[i])*(data[j*N+n]-means[j])*Rinv[i*D+j];
                    }
                }
            #endif  
            clusters->memberships[mem_col_idx+n] = (double)(-0.5f * like + clusters->constant[m] + log(clusters->pi[m])); 
        }
    }
    //finish = clock();
    //DEBUG("estep1: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

void distribution(double* data, clusters_t* clusters, int D, int Mfirst, int Mlast, int N) {
    //clock_t start,finish;
    // Compute likelihood for every data point in each cluster
    double like;
    double* means;
    double* Rinv;
    int D_sq = D*D;
    int R_chn_idx = Mfirst*D_sq;
    int mem_col_idx = Mfirst*N;
    int mea_col_idx = Mfirst*D;
    //start = clock();
    for(int m=Mfirst; m < Mlast; m++,R_chn_idx+=D_sq,mem_col_idx+=N,mea_col_idx+=D) {
        means = (double*) &(clusters->means[mea_col_idx]);
        Rinv = (double*) &(clusters->Rinv[R_chn_idx]);
        double vmax = -1;
        for(int n=0; n < N; n++) {
            like = 0.0;
            #if DIAG_ONLY
                for(int i=0; i < D; i++) {
                    like += (data[i*N+n]-means[i])*(data[i*N+n]-means[i])*Rinv[i*D+i];
                }
            #else
                for(int i=0; i < D; i++) {
                    for(int j=0; j < D; j++) {
                        like += (data[i*N+n]-means[i])*(data[j*N+n]-means[j])*Rinv[i*D+j];
                    }
                }
            #endif  
            double temp =  exp(-0.5 * like );
            if(vmax < temp) vmax = temp;
            clusters->memberships[mem_col_idx+n] = temp; 
        }
        if(vmax<DBL_EPSILON)  vmax = DBL_EPSILON;
        for(int n=0; n < N; n++) {
            clusters->memberships[mem_col_idx+n] /= vmax; 
        }
    }
    //finish = clock();
    //DEBUG("estep1: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

void estep2s(double* data, clusters_t* clusters, int D, int M,int N, double* likelihood) {
    //clock_t start,finish;
    //start = clock();
    double max_likelihood, denominator_sum;
    *likelihood = 0.0;
    for(int n=0; n < N; n++) {
        // initial condition, maximum is the membership in first cluster
        max_likelihood = clusters->memberships[n];
        // find maximum likelihood for this data point
        for(int m=1; m < M; m++) {
            max_likelihood = std::max(max_likelihood,clusters->memberships[m*N+n]);
        }

        // Computes sum of all likelihoods for this event
        denominator_sum = 0.0f;
        for(int m=0; m < M; m++) {
            denominator_sum +=(double) exp(clusters->memberships[m*N+n] - max_likelihood);
        }
        denominator_sum = max_likelihood + (double) log(denominator_sum);
        *likelihood += denominator_sum;

        // Divide by denominator to get each membership
        for(int m=0; m < M; m++) {
            clusters->memberships[m*N+n] = (double) exp(clusters->memberships[m*N+n] - denominator_sum);
            //printf("Membership of event %d in cluster %d: %.3f\n",n,m,clusters->memberships[m*N+n]);
        }
    }
    //finish = clock();
   // DEBUG("estep2: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

void estep2(double* data, clusters_t* clusters, int D, int M,int N, int Nfirst, int Nlast, int tid, double* likelihood) {
    //clock_t start,finish;
    //start = clock();
    double max_likelihood, denominator_sum;
    likelihood[tid] = 0.0;
    for(int n=Nfirst; n < Nlast; n++) {
        // initial condition, maximum is the membership in first cluster
        max_likelihood = clusters->memberships[n];
        // find maximum likelihood for this data point
        for(int m=1; m < M; m++) {
            max_likelihood = std::max(max_likelihood,clusters->memberships[m*N+n]);
        }

        // Computes sum of all likelihoods for this event
        denominator_sum = 0.0f;
        for(int m=0; m < M; m++) {
            denominator_sum +=(double) exp(clusters->memberships[m*N+n] - max_likelihood);
        }
        denominator_sum = max_likelihood + (double) log(denominator_sum);
        likelihood[tid] += denominator_sum;

        // Divide by denominator to get each membership
        for(int m=0; m < M; m++) {
            clusters->memberships[m*N+n] = (double) exp(clusters->memberships[m*N+n] - denominator_sum);
            //printf("Membership of event %d in cluster %d: %.3f\n",n,m,clusters->memberships[m*N+n]);
        }
    }
    //finish = clock();
   // DEBUG("estep2: %f seconds.\n",(double)(finish-start)/(double)CLOCKS_PER_SEC);
}

void mstep_n(double* data, clusters_t* clusters, int D, int Mfirst, int Mlast,  int N) {
    //DEBUG("mstep_n: D: %d, M: %d, N: %d\n",D,M,N);
    int idx_col = Mfirst*N;
    for(int m=Mfirst; m < Mlast; m++,idx_col+=N) {
        clusters->N[m] = 0.0;
        // compute effective size of each cluster by adding up soft membership values
        for(int n=0; n < N; n++) {
            clusters->N[m] += clusters->memberships[idx_col+n];
        }
    }
}

void mstep_mean(double* data, clusters_t* clusters, int D, int Mfirst, int Mlast, int N) {
   // DEBUG("mstep_mean: D: %d, M: %d, N: %d\n",D,M,N);
    int mem_idx_col = Mfirst*N; 
    int mea_idx_col = Mfirst*D;
    for(int m=Mfirst; m < Mlast; m++,mem_idx_col+=N,mea_idx_col+=D) {
        //DEBUG("Cluster %d: ",m);
        int data_idx_col = 0;
        for(int d=0; d < D; d++,data_idx_col+=N) {
            clusters->means[mea_idx_col+d] = 0.0;
            for(int n=0; n < N; n++) {
                clusters->means[mea_idx_col+d] += data[data_idx_col+n]*clusters->memberships[mem_idx_col+n];
            }
            clusters->means[mea_idx_col+d] /= clusters->N[m];
           // DEBUG("%f ",clusters->means[m*D+d]);
        }
        //DEBUG("\n");
    }
}

void mstep_covar(double* data, clusters_t* clusters, int D, int Mfirst, int Mlast, int N) {
    //DEBUG("mstep_covar: D: %d, M: %d, N: %d\n",D,M,N);
    double sum;
    double* means;
    int D_sq = D*D;
    int R_chn_idx = Mfirst*D_sq;
    int mem_col_idx = Mfirst*N;
    for(int m=Mfirst; m < Mlast; m++,R_chn_idx+=D_sq,mem_col_idx+=N) {
        means = &(clusters->means[m*D]);
        int R_row_idx = 0;
        int data_col_idx1 = 0;
        for(int i=0; i < D; i++, R_row_idx+=D,data_col_idx1+=N) {
            int R_col_idx = 0;
            int data_col_idx2 = 0;
            for(int j=0; j <= i; j++,R_col_idx+=D,data_col_idx2+=N) {
                #if DIAG_ONLY
                    if(i != j) {
                        clusters->R[R_chn_idx+R_row_idx+j] = 0.0f;
                        clusters->R[R_chn_idx+R_col_idx+i] = 0.0f;
                        continue;
                    }
                #endif
                sum = 0.0;
                for(int n=0; n < N; n++) {
                    sum += (data[data_col_idx1+n]-means[i])*(data[data_col_idx2+n]-means[j])*clusters->memberships[mem_col_idx+n];
                }
                if(clusters->N[m] > 1.0e-6) {
                    clusters->R[R_chn_idx+R_row_idx+j] = sum / clusters->N[m];
                    clusters->R[R_chn_idx+R_col_idx+i] = sum / clusters->N[m];
                } else {
                    clusters->R[R_chn_idx+R_row_idx+j] = 0.0f;
                    clusters->R[R_chn_idx+R_col_idx+i] = 0.0f;
                }
                if(i==j) clusters->R[R_chn_idx+R_col_idx+i] += 0.001;
            }
        }
    }
}

void* thread_mstep(void* paras)
{
    thread_t* ptr_paras = (thread_t*) paras;

    int first = ptr_paras->td_id * ptr_paras->len;
    int last = first + ptr_paras->len;
    if(last>ptr_paras->M)
        last = ptr_paras->M;
    mstep_n(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    mstep_mean(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    mstep_covar(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    pthread_exit(NULL);
    return NULL;
}

void mt_mstep(double* data, clusters_t* clusters, int D, int M, int N, int t_num)
{
    pthread_t* thread = new pthread_t[t_num];
    pthread_attr_t attr;
    void* status;
    int rc;
    int t;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    thread_t* ptr_tparams = new thread_t[t_num];
    int length =(int) ceil(M / double(t_num));
    for (t = 0; t < t_num; t++) {
        
        ptr_tparams[t].data = data;
        ptr_tparams[t].clusters = clusters;
        ptr_tparams[t].D = D;
        ptr_tparams[t].M = M; ptr_tparams[t].N = N;
        ptr_tparams[t].len = length;
        ptr_tparams[t].td_id =t;
        //printf("Creating thread #%d!\n", t);
        rc = pthread_create(&thread[t], &attr, thread_mstep, (void *)(&ptr_tparams[t]));
        if (rc) {
            //printf("Thread %d: ERROR; return code from pthread_create() is %d\n", t,rc);
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_create");
            //exit(-1);
        }

    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    for (t = 0; t<t_num; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_join");
        }
    }
    delete[] ptr_tparams;
    delete[] thread;
}


void* thread_estep1(void* paras)
{
    thread_t* ptr_paras = (thread_t*) paras;

    int first = ptr_paras->td_id * ptr_paras->len;
    int last = first + ptr_paras->len;
    if(last>ptr_paras->M)
        last = ptr_paras->M;
    estep1(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_n(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_mean(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_covar(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    pthread_exit(NULL);
    return NULL;
}

void mt_estep1(double* data, clusters_t* clusters, int D, int M, int N, int t_num)
{
    pthread_t* thread = new pthread_t[t_num];
    pthread_attr_t attr;
    void* status;
    int rc;
    int t;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    thread_t* ptr_tparams = new thread_t[t_num];
    int length =(int) ceil(M / double(t_num));
    for (t = 0; t < t_num; t++) {
        
        ptr_tparams[t].data = data;
        ptr_tparams[t].clusters = clusters;
        ptr_tparams[t].D = D;
        ptr_tparams[t].M = M; ptr_tparams[t].N = N;
        ptr_tparams[t].len = length;
        ptr_tparams[t].td_id =t;
        //printf("Creating thread #%d!\n", t);
        rc = pthread_create(&thread[t], &attr, thread_estep1, (void *)(&ptr_tparams[t]));
        if (rc) {
            //printf("Thread %d: ERROR; return code from pthread_create() is %d\n", t,rc);
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_create");
            //exit(-1);
        }

    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    for (t = 0; t<t_num; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_join");
        }
    }
    delete[] ptr_tparams;
    delete[] thread;
}

void* thread_estep2(void* paras)
{
    thread_t* ptr_paras = (thread_t*) paras;

    int first = ptr_paras->td_id * ptr_paras->len;
    int last = first + ptr_paras->len;
    if(last>ptr_paras->N)
        last = ptr_paras->N;
    estep2(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, ptr_paras->M, ptr_paras->N, first, last, ptr_paras->td_id, ptr_paras->likelihood);
    //estep2(double* data, clusters_t* clusters, int D, int M, int Nfirst, int Nlast, int tid, double* likelihood)
    //mstep_n(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_mean(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_covar(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    pthread_exit(NULL);
    return NULL;
}

void mt_estep2(double* data, clusters_t* clusters, double* likelihood, int D, int M, int N, int t_num)
{
    pthread_t* thread = new pthread_t[t_num];
    pthread_attr_t attr;
    void* status;
    int rc;
    int t;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    thread_t* ptr_tparams = new thread_t[t_num];
    double* likelihoods = new double[t_num];
    int length =(int) ceil(N / double(t_num));
    for (t = 0; t < t_num; t++) {
        
        ptr_tparams[t].data = data;
        ptr_tparams[t].clusters = clusters;
        ptr_tparams[t].D = D;
        ptr_tparams[t].M = M; ptr_tparams[t].N = N;
        ptr_tparams[t].len = length;
        ptr_tparams[t].td_id =t;
        ptr_tparams[t].likelihood = likelihoods;
        //printf("Creating thread #%d!\n", t);
        rc = pthread_create(&thread[t], &attr, thread_estep2, (void *)(&ptr_tparams[t]));
        if (rc) {
            //printf("Thread %d: ERROR; return code from pthread_create() is %d\n", t,rc);
            delete[] likelihoods;
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_create");
            //exit(-1);
        }

    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    for (t = 0; t<t_num; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            delete[] likelihoods;
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_join");
        }
    }
    *likelihood = 0;
    for(t=0; t<t_num; t++)
        *likelihood += likelihoods[t];
    delete[] likelihoods;
    delete[] ptr_tparams;
    delete[] thread;
}

void* thread_distribution(void* paras)
{
    thread_t* ptr_paras = (thread_t*) paras;

    int first = ptr_paras->td_id * ptr_paras->len;
    int last = first + ptr_paras->len;
    if(last>ptr_paras->M)
        last = ptr_paras->M;
    distribution(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_n(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_mean(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    //mstep_covar(ptr_paras->data, ptr_paras->clusters, ptr_paras->D, first, last,  ptr_paras->N);
    pthread_exit(NULL);
    return NULL;
}

void mt_distribution(double* data, clusters_t* clusters, int D, int M, int N, int t_num)
{
    pthread_t* thread = new pthread_t[t_num];
    pthread_attr_t attr;
    void* status;
    int rc;
    int t;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    thread_t* ptr_tparams = new thread_t[t_num];
    int length =(int) ceil(M / double(t_num));
    for (t = 0; t < t_num; t++) {
        
        ptr_tparams[t].data = data;
        ptr_tparams[t].clusters = clusters;
        ptr_tparams[t].D = D;
        ptr_tparams[t].M = M; ptr_tparams[t].N = N;
        ptr_tparams[t].len = length;
        ptr_tparams[t].td_id =t;
        //printf("Creating thread #%d!\n", t);
        rc = pthread_create(&thread[t], &attr, thread_distribution, (void *)(&ptr_tparams[t]));
        if (rc) {
            //printf("Thread %d: ERROR; return code from pthread_create() is %d\n", t,rc);
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_create");
            //exit(-1);
        }

    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    for (t = 0; t<t_num; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            delete[] ptr_tparams;
            delete[] thread;
            mexErrMsgTxt("ERROR return code from pthread_join");
        }
    }
    delete[] ptr_tparams;
    delete[] thread;
}



#endif

