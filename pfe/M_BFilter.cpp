/*************************************************************
*    Piecewise flat embedding for image segmentation
*    Bilateral Filter
*
*    Author: Chaowei FANG
*    Date: 12/3/2015
*
*******************************************************************/
 
#include "mex.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )     
{    
    int v,i,j,k,ni,nj,ndim,rows,cols;
	int px,channels;
	int ii,jj;
    double *ima,*fima,*weight,sigI,sigS,w;
    const int  *dims ;
    
    if(nrhs<4)
    {
        mexErrMsgTxt("Error: Wrong number of arguments!! must be 4\n");
        //return;
    } 
   
    /*Get image*/
    ima = (double*)mxGetData(prhs[0]);    
    //rows=(int)mxGetM(prhs[0]);
    //cols=(int)mxGetN(prhs[0]);
	
    ndim = (int)mxGetNumberOfDimensions(prhs[0]);
	
    dims = mxGetDimensions(prhs[0]);
	
	if(ndim != 2 && ndim!=3)
	{
		mexErrMsgTxt("The dimension of input image is gray or color");
	}
    px = dims[0]*dims[1];
	rows = dims[0];
	cols = dims[1];
	if(ndim == 2)
		channels = 1;
	else
		channels = (int)dims[2];
	//printf("dims=%d %d %d",dims[0],dims[1],dims[2]);
    /*Get the Integer*/
    v = (int)(mxGetScalar(prhs[1]));
    sigS = mxGetScalar(prhs[2]);  
    sigI = mxGetScalar(prhs[3]);  
    
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
    fima = (double*)mxGetData(plhs[0]);
	weight = (double*)mxMalloc(px*sizeof(double));
    //weight = mxGetPr( mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL));
    printf("rows=%d cols=%d\n",rows,cols);
	//Set initialization values
    for(i=0;i<px;i++)
    {

		for(j=0; j<channels; j++)
			fima[px*j+i] = ima[px*j+i];
        weight[i]=1;
    }

	//compute bilateral filter
	for(i=0;i<rows;i++) 
    for(j=0;j<cols;j++) 
    {
        for(ii=0;ii<=v;ii++) 
        for(jj=-v;jj<=v;jj++)
        {
            ni=i+ii;
            nj=j+jj;
            if(ii==0 && jj==0) continue;
            if(ni>=0 && nj>=0 && ni<rows && nj<cols)
            {      
				d = 0;
				for(k=0;k<channels;k++)
                {
					double tempd=(ima[nj*rows+ni+k*px]-ima[j*rows+i+k*px]);
					d +=tempd*tempd;
				}
                w=exp(-(ii*ii+jj*jj)/(2*sigS*sigS)+(-d/(2*sigI*sigI)));
                //if(i==0 && j==0)
                //    printf("sigS=%f sigI=%f d=%f w=%f\n",sigS,sigI,d,w);
                // symmetric weight computation
                for(k=0;k<channels;k++)
				{
					fima[j*rows+i + k*px]  += w*ima[nj*rows+ni+k*px];
					fima[nj*rows+ni+k*px]  += w*ima[j*rows+i+k*px];
				}

                weight[j*rows+i]+= w;
				weight[nj*rows+ni]+= w;
            }
        }
    }

    for(i=0;i<px;i++) 
	{
		for(j=0; j<channels; j++)
		   fima[i+j*px]/=weight[i];   
	}
	mxFree(weight);
    //return;
}
