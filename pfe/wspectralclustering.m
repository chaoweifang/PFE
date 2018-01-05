function y = wspectralclustering(I,J,val,nls,Chs)

W=sparse(I,J,val);
S = full(sum(W, 1));
[wx, wy] = size(W);
x = 1 : wx;
D = sparse(x, x, S, wx, wy);
opts.issym=1;
opts.isreal = 1;
opts.disp=0;
[EigVectLE, EigVal] = eigs(D - W, D, Chs+1, 'sm',opts);
EigVectLE = EigVectLE(:,1:Chs);  
EigVal = diag(EigVal)';
EigVal = EigVal(1:Chs);
flagp = EigVal>0;
EigVal = EigVal(flagp);
EigVect = EigVectLE(:,flagp);
EigVect = bsxfun(@rdivide,EigVect,sqrt(EigVal));
[~,y]=vl_kmeans(EigVect',nls,'MaxNumIterations',500,'NumRepetitions',10);
y = y';

end