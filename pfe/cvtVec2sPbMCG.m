function sPb = cvtVec2sPbMCG(vec,rs,cs)
    vminmax = minmax(vec');
    vmin = vminmax(:,1); vmax=vminmax(:,2);
    vmm = vmax-vmin;
    zero_vmm = vmm==0;
    vmm(zero_vmm)=1;
    vec = bsxfun(@rdivide,bsxfun(@minus,vec,vmin'),vmm');
    %if vmm==0, vmm=1; end
    %vec = (vec-vmin)/vmm;
    
    nvec = size(vec,2);
    vect = reshape(vec,rs,cs,nvec);
    EigVal = ones(1,nvec)/nvec;

    sPb_thin = zeros(2*rs+1, 2*cs+1);
    for v = 1 : nvec
            vectemp = vect(:,:,v)/sqrt(EigVal(v));
            sPb_thin = sPb_thin + seg2bdry_wt(vectemp, 'doubleSize');
         %end
    end
    sPb = sPb_thin.^(1/sqrt(2));
end