%Piecewise flat embedding for image segmentation
%Erase small regions with tiny width. 
%@Chaowei FANG
%12/3/2015

function result = M_BFilter(image, radius, sigmaS,sigmaI)
%in
%image     m*n*d double   image(m*n pixels) data
%radius    double         filter window size
%sigmaS    double         sigma of spatial filtering
%sigmaI    double         sigma of intensity filtering

%out
%result    m*n*d          double filtering results
