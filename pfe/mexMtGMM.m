%Piecewise flat embedding for image segmentation
%GMM Estimation Using multithreads
%@Chaowei FANG
%12/3/2015

%function [clusters,covRs,phi,probs,logl] = mexMtGMM(data,nclusters, nthreads, inicenters,inicovRs)
%clusters number: m 
%data number: n, dimension: d
%in
%data d*n        input data
%nclusters m     number of clusters
%inicenters d*m  cluster centers initialization
%iicovRs   d*d*m covariance matrix initialization

%out
%clusters d*m  cluster centers
%covRs d*d*m   covariance matrix
%phi 1*m       weights of clusters
%probs n*m     probabilites to clusters of data
%logl  1       log of likelihood of GMM
