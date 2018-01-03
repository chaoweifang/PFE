function Yfinal= NestedBregmanPRW(I,J,Wval, Yin, p, miu, r,rex,rscale,rlb,sigma,tol)
if nargin < 8
    disp('8 arguments are required at least!');
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = eps;
end
if ~exist('tol','var') || isempty(tol)
    tol = eps;
end
if ~exist('rscale','var') || isempty(rscale)
    rscale = 100;
end
if ~exist('rlb','var') || isempty(rlb)
    rlb = 1e-6;
end
%
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
clear vidx_M

[total,m] = size(Yin);
W = sparse(I,J,Wval,total,total);
clear I J Wval

d = full(sum(W,2));
clear W

dsq = sqrt(d);
D   = sparse(1:total,1:total,d);
S = miu*Mt*M;
SD = S+r*D;
[L_SD,~,P_SD]=lchol(SD);
clear SD

Orthog = bsxfun(@times,Yin,dsq);
[tempu,~,tempv] = svd(Orthog,0);
P0 = tempu*tempv';
clear Orthog tempu tempv

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

    %disp(['inner it:' num2str(count) ' v=' num2str(v1)])
    T1 = bsxfun(@times,F,r*dsq);
    while err_l1 >tol && count<8
        T=T1+Mt*(miu*(DiffY0-DiffBY0));
        T = T(P_SD,:);
        Y = lsolve(L_SD,T);
        Y(P_SD,:) = Y;
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

%disp('stage 2');
clear L_SD P_SD SD u v DiffY0 DiffBY0 T1 T Y
SD = S+rex*D;
clear S

[L_SD,~,P_SD]=lchol(SD);
clear SD

iters=0;
M0 = M; Mt0=Mt; Widx0 = Widx;


err_lp =tol+100;
F=P0-DiffBP0;
clear P0 DiffBP0

T1=bsxfun(@times,F,rex*dsq);
v_np1 = rex/2*sum(sum(F.^2));
%disp(['outer it:' num2str(iters) ' v=' num2str(v_np1)])

while err_lp>tol && iters<5   
    DiffY0 = zeros(prnum,m);
    DiffBY0 = zeros(prnum,m);
    v1 = rex/2*sum(sum(F.^2));
    %disp(['inner it:' num2str(count) ' v=' num2str(v1)]); 
    err_l1 = 2*tol;
    count=0;
    while err_l1 >tol && count<20
        T=T1+Mt*(miu*(DiffY0-DiffBY0));
        T = T(P_SD,:);
        Y = lsolve(L_SD,T);
        Y(P_SD,:) = Y;
        [DiffY0,DiffBY0,v2] = mexMtShrink(Iidx,Jidx,Widx,Y,DiffBY0,1/miu,8);
        v_or = rex/2* sum(sum((bsxfun(@times,Y,dsq)-F).^2));
        v2 = v2 + v_or;
        vtemp =v1+v2;
        err_l1 = abs((v1-v2)/(vtemp+(vtemp==0)));    
        v1 = v2;
        count=count+1;
        %disp(['inner it:' num2str(count) ' err_l1 = ' num2str(err_l1) ' v=' num2str(v2)])
        %disp(['err_l1 = ' num2str(err_l1)])
    end
    clear DiffY0 DiffBY0 T L_SD P_SD
    
    %energy calculation
    diffY =abs(M0*Y);
    res = sum(diffY,2);
    v_np2 = sum(res.^p)+v_or;
    v_temp = v_np1+v_np2;
    err_lp = abs(v_np1-v_np2)/(v_temp+(v_temp==0));
    v_np1 = v_np2;
    %disp(['outer it:' num2str(iters) ' err_l1 = ' num2str(err_lp) ' v=' num2str(v_np2)])
    %disp('-------------------------------------------------------')

    %normalize
    vmin = min(Y,[],1);
    vmax = max(Y,[],1);
    vmm = vmax-vmin;
    vmm = vmm+(vmm==0);
    Y_norm = bsxfun(@rdivide,bsxfun(@minus,Y,vmin),vmm);

    %reweight
    tmp = sqrt(sum(bsxfun(@times,Y,dsq).^2,1)+eps);
    EigVal = sum(diffY,1)./tmp;
    EigVal = EigVal + (EigVal==0);
    Y = bsxfun(@rdivide,Y_norm, sqrt(EigVal));
    clear Y_norm

    iters =iters+1;
    if err_lp<=tol || iters>=5
        break;
    end
    if p~=1
        res = sum(abs(M0*Y),2);
        r_p = max(rscale*res,rlb).^(p-1);
        Widx = Widx0.*r_p;
        R = sparse(1:prnum,1:prnum,r_p);
        clear r_p Yt;
        Mt = Mt0*R;
        M = R*M0;
        SD = miu*Mt*M+rex*D;
        [L_SD,~,P_SD]=lchol(SD);
        clear SD R
    else
        break;
    end 
end

Yfinal = Y;

end
