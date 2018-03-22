#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "invert_matrix.h"
#include "gaussian.h"


static double double_abs(double x);

static int 
ludcmp(double *a,int n,int *indx,double *d);

static void 
lubksb(double *a,int n,int *indx,double *b);

/*
 * Inverts a square matrix (stored as a 1D double array)
 * 
 * actualsize - the dimension of the matrix
 *
 * written by Mike Dinolfo 12/98
 * version 1.0
 */
void invert_cpu(double* data, int actualsize, double* log_determinant)  {
    int maxsize = actualsize;
    int n = actualsize;
    *log_determinant = 0.0;
    
    if (actualsize == 1) { // special case, dimensionality == 1
        *log_determinant = logf(data[0]);
        data[0] = 1.0f / data[0];
    } else if(actualsize >= 2) { // dimensionality >= 2
        for (int i=1; i < actualsize; i++) data[i] /= data[0]; // normalize row 0
        for (int i=1; i < actualsize; i++)  { 
            for (int j=i; j < actualsize; j++)  { // do a column of L
                double sum = 0.0;
                for (int k = 0; k < i; k++)  
                    sum += data[j*maxsize+k] * data[k*maxsize+i];
                data[j*maxsize+i] -= sum;
            }
            if (i == actualsize-1) continue;
            for (int j=i+1; j < actualsize; j++)  {  // do a row of U
                double sum = 0.0;
                for (int k = 0; k < i; k++)
                    sum += data[i*maxsize+k]*data[k*maxsize+j];
                data[i*maxsize+j] = 
                    (data[i*maxsize+j]-sum) / data[i*maxsize+i];
            }
        }

        for(int i=0; i<actualsize; i++) {
            *log_determinant += (double) log(fabs(data[i*n+i]));
        }
        //printf("\n\n");
        for ( int i = 0; i < actualsize; i++ )  // invert L
            for ( int j = i; j < actualsize; j++ )  {
                double x = 1.0;
                if ( i != j ) {
                    x = 0.0;
                    for ( int k = i; k < j; k++ ) 
                        x -= data[j*maxsize+k]*data[k*maxsize+i];
                }
                data[j*maxsize+i] = x / data[j*maxsize+j];
            }
        for ( int i = 0; i < actualsize; i++ )   // invert U
            for ( int j = i; j < actualsize; j++ )  {
                if ( i == j ) continue;
                double sum = 0.0;
                for ( int k = i; k < j; k++ )
                    sum += data[k*maxsize+j]*( (i==k) ? 1.0f : data[i*maxsize+k] );
                data[i*maxsize+j] = -sum;
            }
        for ( int i = 0; i < actualsize; i++ )   // final inversion
            for ( int j = 0; j < actualsize; j++ )  {
                double sum = 0.0;
                for ( int k = ((i>j)?i:j); k < actualsize; k++ )  
                    sum += ((j==k)?1.0f:data[j*maxsize+k])*data[k*maxsize+i];
                data[j*maxsize+i] = sum;
            }
    } else {
        //PRINT("Error: Invalid dimensionality for invert(...)\n");
    }
 }


