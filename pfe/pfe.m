% ------------------------------------------------------------------------
function vecs = pfe(lbp,Yor,param_s,nvec)

[tx,ty]=size(lbp);

% build the pairwise affinity matrix        
l{1} = zeros(tx + 1, ty);
l{1}(2:end, :) = lbp;
l{2} = zeros(tx, ty + 1);
l{2}(:, 2:end) = lbp;
[val,I,J] = buildW(l{1},l{2}, param_s.dthresh, param_s.ic_gamma);
%piecewise flat embedding
if param_s.p~=1
    if strcmp(param_s.wtype,'w')
        if strcmp(param_s.version,'iccv')
            EigVect = NestedBregmanPRW1(val,I,J,Yor,param_s.p,param_s.miu,param_s.r,param_s.rex,param_s.rlb);
        elseif strcmp(param_s.version,'journal')
            EigVect = NestedBregmanPRW(val,I,J,Yor,param_s.p,param_s.miu,param_s.r,param_s.rex,param_s.rlb);
        else
            error(['wrong version ' param_s.version])
        end
        EigVal = [];
    elseif strcmp(param_s.wtype,'o')
        [EigVect,EigVal] = NestedBregmanP(val,I,J,Yor,param_s.p,param_s.miu,param_s.r,param_s.rex,param_s.rlb);
    else
        error(['error weighting type ' param_s.wtype])
    end
else
    if strcmp(param_s.version,'iccv')
        [EigVect,EigVal]= NestedBregman1(val,I,J, Yor, param_s.miu, param_s.r,param_s.rex);
    elseif strcmp(param_s.version, 'journal')
        [EigVect,EigVal]= NestedBregman(val,I,J, Yor, param_s.miu, param_s.r,param_s.rex);
    else
        error(['wrong version ' param_s.version])
    end
end
clear val I J Yor;
for v = 1 : nvec
    temp = reshape(EigVect(:, v), [ty tx])';
    EigVect(:,v) = temp(:);
end
vecs.vect=EigVect;
vecs.val =EigVal;
end
