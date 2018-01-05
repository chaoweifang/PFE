function sPb = cvtVec2Edge(vec,im_r,im_c)

    nvec = size(vec,2);
    vect = vec;
    vect = reshape(vect,im_r,im_c,nvec);
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
            vec = vect(:,:,v);%/sqrt(EigVal(v));
            for o = 1 : norient
                theta = dtheta*o;
                f = oeFilter(sigma, support, theta, deriv, hil);
                sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, vec));
            end
    end
    sPb = max(sPb,[],3).^(1/sqrt(2));
end
