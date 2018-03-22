#ifndef INVERT_MATRIX_H
#define INVERT_MATRIX_H

void invert_cpu(double* data, int actualsize, double* log_determinant);
int invert_matrix(double* a, int n, double* determinant);
#endif

