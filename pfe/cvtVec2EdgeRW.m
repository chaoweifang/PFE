function sPb = cvtVec2EdgeRW(vecs,im_r,im_c)

    
    vect = vecs.vect;
    nvec = size(vect,2);
    vmin = min(vect,[],1);
    vmax = max(vect,[],1);
    vmm = vmax-vmin;
    vmm = vmm+(vmm==0);
    vect = bsxfun(@rdivide,bsxfun(@minus,vect,vmin),vmm);
    vect = reshape(vect,im_r,im_c,nvec);
    %vminmax= minmax(vect(:)');
    %vmin = vminmax(1);vmax = vminmax(2);
    %vmm = vmax-vmin;
    %vect = (vect-vmin)/(vmm+(vmm==0));
    %vect = reshape(vect,im_r,im_c,nvec);
    % OE parameters
    hil = 0;
    deriv = 1;
    support = 3;
    sigma = 1;
    norient = 8;
    dtheta = pi/norient;
    ch_per = [4 3 2 1 8 7 6 5];

    sPb = zeros(im_r, im_c, norient);
    for v = 1 : nvec
        %if EigVal(v) > 0,
            vec = vect(:,:,v)*vecs.val(v);
            for o = 1 : norient,
                theta = dtheta*o;
                f = oeFilter(sigma, support, theta, deriv, hil);
                sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, vec));
            end
        %end
    end
    
    sPb = max(sPb,[],3).^(1/sqrt(2));
    
end