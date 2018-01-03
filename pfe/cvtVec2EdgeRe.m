function sPb = cvtVec2EdgeRe(vecs,im_r,im_c)
%     if nargin < 4 || isempty(val)
%         EigVal = ones(1,size(vec,2))/size(vec,2);
%     else
%         EigVal = val;
%     end

    
    vect = vecs.vect;
    nvec = size(vect,2);
    %vminmax= minmax(vect(:)');
    %vmin = vminmax(1);vmax = vminmax(2);
    %vmm = vmax-vmin;
    %vect = (vect-vmin)/(vmm+(vmm==0));
    vmin = min(vect,[],1);
    vmax = max(vect,[],1);
    vmm = vmax-vmin;
    vmm = vmm+(vmm==0);
    vect = bsxfun(@rdivide,bsxfun(@minus,vect,vmin),vmm);
    vect = reshape(vect,im_r,im_c,nvec);
    % OE parameters
    hil = 0;
    deriv = 1;
    support = 3;
    sigma = 1;
    norient = 8;
    dtheta = pi/norient;
    ch_per = [4 3 2 1 8 7 6 5];
    
    vals = vecs.val;
    sPb = zeros(im_r, im_c, norient);
    %wrs = 1./(EigVal+(EigVal==0));
    %wrs = wrs/sum(wrs);
    for v = 1 : nvec
        %if EigVal(v) > 0,
            vec = vect(:,:,v)/sqrt(vals(v)+(vals(v)==0));
            for o = 1 : norient,
                theta = dtheta*o;
                f = oeFilter(sigma, support, theta, deriv, hil);
                sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, vec));
            end
        %end
    end
    
    sPb = max(sPb,[],3).^(1/sqrt(2));
    
end
