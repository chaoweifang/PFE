function [Yfinal,EigVal]= NestedBregman1(I,J,Wval, Yin,miu, r,rex,sigma,tol)
if nargin < 7
    disp('7 arguments are required at least!');
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = eps;
end
if ~exist('tol','var') || isempty(tol)
    tol = eps;
end

pridx = J>I;
Iidx = I(pridx);
Jidx = J(pridx);
Widx = Wval(pridx);
prnum = numel(Iidx);
ridx_M = repmat((1:prnum)',2,1);
cidx_M = [Iidx;Jidx];
vidx_M = [Widx;-Widx];
M = sparse(ridx_M(:),cidx_M(:),vidx_M(:));
Mt= sparse(cidx_M(:),ridx_M(:),vidx_M(:));
[total,m] = size(Yin);
W = sparse(I,J,Wval,total,total);
d = full(sum(W,2));
dsq = sqrt(d);
D   = sparse(1:total,1:total,d);
S = miu*Mt*M;
SD = S+r*D;
%[LD_SD,~,P_SD]=ldlchol(SD);

Orthog = bsxfun(@times,Yin,dsq);
[tempu,~,tempv] = svd(Orthog,0);
P0 = tempu*tempv';
Y = bsxfun(@rdivide, P0, dsq);
DiffBP0 = zeros(total,m);%error back of p   

k=0;
err_ort = 2*sigma;
while err_ort>sigma && k<5
    F = P0-DiffBP0;
    err_l1 = 2*tol;
    count=0;
    DiffY0 = zeros(prnum,m);
    DiffBY0 = zeros(prnum,m);
    v1 = r/2*sum(sum(F.^2));
    %v1 = r/2*sum(sum((bsxfun(@times,Y,dsq)-F).^2));
    %v1 = sum(sum(abs(diffY))) + r/2*sum(sum((bsxfun(@times,Y,dsq)-F).^2));
    %disp(['inner it:' num2str(count) ' v=' num2str(v1)])
    T1 = bsxfun(@times,F,r*dsq);
    while err_l1 >tol && count<8
        T=T1+Mt*(miu*(DiffY0-DiffBY0));
        Y=SD\T;
        %T = T(P_SD,:);
        %Y = ldlsolve(LD_SD,T);
        %Y(P_SD,:) = Y;
        %diffY  = M*Y;
        %v2 = sum(sum(abs(diffY))) + r/2* sum(sum((bsxfun(@times,Y,dsq)-F).^2));
        [DiffY0,DiffBY0,v2] = mexMtShrink(Iidx,Jidx,Widx,Y,DiffBY0,1/miu,8);
        v2 = v2 + r/2* sum(sum((bsxfun(@times,Y,dsq)-F).^2));
        vtemp =v1+v2;
        err_l1 = abs((v1-v2)/(vtemp+(vtemp==0))); 
        v1 = v2;
        count=count+1;
        %disp(['inner it:' num2str(count) ' err_l1 = ' num2str(err_l1) ' v=' num2str(v2)])   
    end
    X = bsxfun(@times,Y,dsq) + DiffBP0;
    [u,~,v] = svd(X,0);
    P0 = u*v';
    DiffBP0 = X-P0;
    err_ort = sum(sum((Y'*D*Y-eye(m)).^2));
    %disp(['err_ort = ' num2str(err_ort)])
    %disp('-----------------------')
    k=k+1;
end
%Ystep1 = Y;

%clear LD_SD P_SD SD
SD = S+rex*D;

%[LD_SD,~,P_SD]=ldlchol(SD);
err_l1 = 2*tol;
count=0;
DiffY0 = zeros(prnum,m);
DiffBY0 = zeros(prnum,m);

F=P0-DiffBP0;
T1=bsxfun(@times,F,rex*dsq);
v1 = rex/2*sum(sum(F.^2));
%disp(['inner it:' num2str(count) ' v=' num2str(v1)])
%v1 = rex/2*sum(sum((Y-Ystep1).^2));
while err_l1 >tol && count<40
    T=T1+Mt*(miu*(DiffY0-DiffBY0));
    %T = T(P_SD,:);
    %Y = ldlsolve(LD_SD,T);
    %Y(P_SD,:) = Y;
    Y = SD\T;
    %diffY  = M*Y;
    %v2 = sum(sum(abs(diffY))) + rex/2* sum(sum((bsxfun(@times,Y,dsq)-F).^2));
    [DiffY0,DiffBY0,v2] = mexMtShrink(Iidx,Jidx,Widx,Y,DiffBY0,1/miu,8);
    v2 = v2 + rex/2* sum(sum((bsxfun(@times,Y,dsq)-F).^2));
    vtemp =v1+v2;
    err_l1 = abs((v1-v2)/(vtemp+(vtemp==0)));    
    v1 = v2;
    count=count+1;
    %disp(['inner it:' num2str(count) ' err_l1 = ' num2str(err_l1) ' v=' num2str(v2)])
    %disp(['err_l1 = ' num2str(err_l1)])

end
Yfinal = Y;
tmp = sqrt(sum((bsxfun(@times,Yfinal,dsq)).^2,1)+eps);
EigVal = sum(abs(M*Yfinal),1)./tmp;
end
