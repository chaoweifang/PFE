// IndexSub.cpp : Defines the entry point for the console application.
//

#include "mex.h"
#include "math.h"
#include "pthread.h"
#include "string.h"
//# include <memory.h>
// used for accomplishing substraction operation C=A-B
pthread_mutex_t mutexev;
pthread_mutex_t mutexerr;
//Parameters required to transit into single threads
typedef struct
{
	double shrpara;
	double *affi_r;
	double *affi_c;
	double *affi_val;
	double *Y;
	double *diff_I_BY0;
	double *ev;
	double *diff_Y0;
	double *diff_O_BY0;
	//double* matrix;  //used to transit the matrix Minuend
	//double shrinkfactor;
	//double* matrix2;  //used to transit the matrix Subtrahend
	//double* factors; //used to transit the array multiplier
	//double* data;    //used to transit the mulitiplication result
	ulong_T diff_rs;       //rows of diff matrices and affinity arrays
	ulong_T diff_cs;       //columns of diff matrices
	ulong_T Y_rs;       //rows of Y
	ulong_T Y_cs;       //columns of Y
	uint8_T t_id;     //present thread id (0...n-1)
	ulong_T t_length; //length each thread processed
	bool* errflag;
} ThreadParas, *PThreadParas;

//Single thread to calculate bsxfunction
void *ThreadCal(void *para)
{
	//int i;
	//long tid;
	//double result = 0.0;
	PThreadParas ptr_tparams = (PThreadParas)para;
	double shrpara=ptr_tparams->shrpara;
	double* affi_r=ptr_tparams->affi_r;
	double* affi_c=ptr_tparams->affi_c;
	double*	affi_val=ptr_tparams->affi_val;
	double*	Y=ptr_tparams->Y;
	double*	diff_I_BY0=ptr_tparams->diff_I_BY0; 
	double* ev=ptr_tparams->ev;
	double*	diff_Y0=ptr_tparams->diff_Y0;
	double* diff_O_BY0=ptr_tparams->diff_O_BY0; //output data
	ulong_T	Y_rs=(ulong_T)ptr_tparams->Y_rs,       Y_cs=(ulong_T)ptr_tparams->Y_cs;
	ulong_T	diff_rs=(ulong_T)ptr_tparams->diff_rs, diff_cs=(ulong_T)ptr_tparams->diff_cs;
	bool* errflag = ptr_tparams->errflag;
	ulong_T r_start = ptr_tparams->t_id * ptr_tparams->t_length;
	ulong_T r_end = r_start + ptr_tparams->t_length;
	if (r_end > (ulong_T)diff_rs) r_end = (ulong_T)diff_rs;
	double  val_affi_temp;
	double diff_idx_affi;
	double diff_temp;
	ulong_T idx;
	ulong_T j;
	ulong_T idx_affi_r;
	ulong_T idx_affi_c;
	ulong_T i;
	double ev_temp = 0.0l;
	//mexPrintf("thread #%d: start=%d end=%d\n", ptr_tparams->t_id,r_start,r_end);
	//double* matrix= ptr_tparams->matrix;
	//double shrpara = ptr_tparams->shrinkfactor;
	//double* matrix2= ptr_tparams->matrix2;
	//double* factors= ptr_tparams->factors;
	//double* data = ptr_tparams->data;
	//printf("Starting thread #%d!\n", ptr_tparams->t_id);
	//printf("r_end = %d",r_end);
	for (i = r_start; i < r_end; i++)
	{
		//get indexes of pixel pair in affinity array
		idx_affi_r = (ulong_T)affi_r[i]-1;
		idx_affi_c = (ulong_T)affi_c[i]-1;
		val_affi_temp = affi_val[i];
		//check if idexes are out of boundary
		if(idx_affi_r>=Y_rs || idx_affi_c>=Y_rs)
		{
		    pthread_mutex_lock (&mutexerr);
		    *errflag = true;
		    pthread_mutex_unlock (&mutexerr);
		    pthread_exit(NULL);
		}

		for (j = 0,idx=i; j < diff_cs; j++,idx += diff_rs,idx_affi_c += Y_rs, idx_affi_r += Y_rs)
		{
			//idx = j*diff_rs + i; //get index of diff matrices

		    //affinity weighted difference between 2 pixels
		    diff_idx_affi =  val_affi_temp * (Y[idx_affi_r] - Y[idx_affi_c]);
		    diff_temp = diff_idx_affi + diff_I_BY0[idx];

		    //shinkage
			if(diff_temp <= shrpara && diff_temp >= -shrpara)
				diff_Y0[idx] = 0.0l;
			else if(diff_temp > shrpara)
				diff_Y0[idx] = diff_temp - shrpara;
			else
				diff_Y0[idx] = diff_temp + shrpara;

			diff_O_BY0[idx] = diff_temp - diff_Y0[idx];
			ev_temp +=  diff_idx_affi>0? diff_idx_affi:-diff_idx_affi;
			//if(diff_idx_affi< 0.0l)
			//	ev_temp -=  diff_idx_affi;
			//else
			//	ev_temp +=  diff_idx_affi;

			//calculate energy value
			//pthread_mutex_lock (&mutexev);
			//ptr_tparams->ev[ptr_tparams->t_id] = ptr_tparams->ev[ptr_tparams->t_id]  + abs(diff_idx_affi);
			//if(i==2)
			//	printf("abs(diff_idx_affi)=%f ",abs(diff_idx_affi));
			//pthread_mutex_unlock (&mutexev);
		}
		//printf("Thread #%d: processing line %d!\n",ptr_tparams->t_id,i);
	}
	ptr_tparams->ev[ptr_tparams->t_id] = ev_temp;
	//printf("\n");
	//printf("Thread %ld eval=%f\n", ptr_tparams->t_id,ptr_tparams->ev[ptr_tparams->t_id]);
	/*for (i = 0; i<1000000; i++)
	{
		result = result + sin(i) * tan(i);
	}
	printf("Thread %ld done. Result = %e\n", tid, result);*/
	pthread_exit(NULL);
	return NULL;
}



void mexFunction(
	int nargout,
	mxArray *out[],
	int nargin,
	const mxArray *in[]
	)
{
	/* declare variables */
	//size_t nc, nv;
	//size_t j;
	//double *A, *C; 
	int t_num = 0;   //nubmer of threads
	double shrpara, *affi_r, *affi_c, *affi_val, *Y, *diff_I_BY0; //input data
	double *ev, *diff_Y0, *diff_O_BY0; //output data
	bool errflag=false;

	//ulong_T total;
	ulong_T diff_rs, diff_cs, Y_rs, Y_cs;

	if (nargin < 6 || nargin>7)
		mexErrMsgTxt("6 or 7 parameters are required!");

	if(!mxIsDouble(in[0]) ||  !mxIsDouble(in[1]) || !mxIsDouble(in[2]) || !mxIsDouble(in[3]) || !mxIsDouble(in[4]))
		mexErrMsgTxt("All input data required to be double!");

	//get size of diff array
	diff_rs = (ulong_T)(mxGetM(in[4]));
	diff_cs = (ulong_T)(mxGetN(in[4]));

	if(diff_rs != (ulong_T)mxGetM(in[0]) || diff_rs != (ulong_T)mxGetM(in[1]) || diff_rs != (ulong_T)mxGetM(in[2]))
		mexErrMsgTxt("Input affity arrays required in same rows as diff matrix");

	//get size of Y
	Y_rs = (ulong_T)(mxGetM(in[3]));
	Y_cs = (ulong_T)(mxGetN(in[3]));
	if(Y_cs != diff_cs)
		mexErrMsgTxt("Rows of Y and diff matrix must be equal.");
	//Get input data required for shrinking
    affi_r   = mxGetPr(in[0]);
    affi_c   = mxGetPr(in[1]);
    affi_val = mxGetPr(in[2]);
    Y        = mxGetPr(in[3]);
    diff_I_BY0 = mxGetPr(in[4]);
	shrpara  = mxGetScalar(in[5]);

	//get number of threads
	if (nargin <7)
		t_num = 1;
	else
		t_num = (int)mxGetScalar(in[6]);

	//Create the output matrix
	out[0] = mxCreateNumericMatrix(diff_rs, diff_cs, mxDOUBLE_CLASS, mxREAL);
	diff_Y0 = mxGetPr(out[0]);
	out[1] = mxCreateNumericMatrix(diff_rs, diff_cs, mxDOUBLE_CLASS, mxREAL);
	diff_O_BY0 = mxGetPr(out[1]);
	out[2] = mxCreateDoubleScalar(0);
	ev = mxGetPr(out[2]);

	/* Initialize and set thread detached attribute */
	pthread_t* thread = new pthread_t[t_num];
	double* ev_ths = new double[t_num];
	memset(ev_ths,0,t_num*sizeof(double));
	pthread_attr_t attr;
	void* status;
	int rc;
	int t;
	//pthread_mutex_init(&mutexev, NULL);
	pthread_mutex_init(&mutexerr, NULL);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
	//Create n threads
	PThreadParas ptr_tparams = new ThreadParas[t_num];
	ulong_T length =(ulong_T) ceil(diff_rs / double(t_num));
	for (t = 0; t < t_num; t++) {
		//input data
		ptr_tparams[t].shrpara = shrpara;
		ptr_tparams[t].affi_r = affi_r;
		ptr_tparams[t].affi_c = affi_c;
		ptr_tparams[t].affi_val = affi_val;
		ptr_tparams[t].Y = Y;
		ptr_tparams[t].diff_I_BY0 = diff_I_BY0; 
		ptr_tparams[t].ev = ev_ths;
		ptr_tparams[t].diff_Y0 = diff_Y0;
		ptr_tparams[t].diff_O_BY0 = diff_O_BY0; //output data
		//ptr_tparams[t].matrix = A;
		//ptr_tparams[t].shrinkfactor = shrpara;
		//ptr_tparams[t].data = C;
		ptr_tparams[t].Y_rs = Y_rs;       ptr_tparams[t].Y_cs = Y_cs;
		ptr_tparams[t].diff_rs = diff_rs; ptr_tparams[t].diff_cs = diff_cs;
		ptr_tparams[t].t_length = length; ptr_tparams[t].t_id =(uint8_T)t;
		ptr_tparams[t].errflag = &errflag;
		//mexPrintf("Creating thread #%d!\n", t);
		rc = pthread_create(&thread[t], &attr, ThreadCal, (void *)(&ptr_tparams[t]));
		if (rc) {
			//printf("Thread %d: ERROR; return code from pthread_create() is %d\n", t,rc);
			//pthread_mutex_destroy(&mutexev);
			pthread_mutex_destroy(&mutexerr);
			delete[] ptr_tparams;
			delete[] thread;
			delete[] ev_ths;
			mexErrMsgTxt("ERROR return code from pthread_create");
			//exit(-1);
		}

	}

	/* Free attribute and wait for the other threads */
	pthread_attr_destroy(&attr);
	for (t = 0; t<t_num; t++) {
		rc = pthread_join(thread[t], &status);
		if (rc) {
			pthread_mutex_destroy(&mutexerr);
			delete[] ptr_tparams;
			delete[] thread;
			delete[] ev_ths;

			mexErrMsgTxt("ERROR return code from pthread_join");
		}
	}
	pthread_mutex_destroy(&mutexerr);
	delete[] ptr_tparams;
	delete[] thread;

	for(t=0; t<t_num; t++)
	{
		*ev = *ev + ev_ths[t];
		//mexPrintf("ev = %f\n",ev_ths[t]);
	}
	delete[] ev_ths;
	//pthread_mutex_destroy(&mutexev);
			
	if(errflag)
		mexErrMsgTxt("Index array of affinity out of boundary");
}