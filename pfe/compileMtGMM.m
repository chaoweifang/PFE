%mex -largeArrayDims mexMtGMM.cpp invert_matrix.cpp -I'E:\toolbox\pthreadswin32\Pre_built\include' -L'E:\\toolbox\\pthreadswin32\\Pre_built\\lib\\x64' -lpthreadVC2
mex -largeArrayDims mexMtGMM.cpp invert_matrix.cpp -lpthread
mex -largeArrayDims mexGaussianPdf.cpp invert_matrix.cpp -lpthread
mex -largeArrayDims mexGaussianPdf1.cpp invert_matrix.cpp -lpthread -outdir ../lib/