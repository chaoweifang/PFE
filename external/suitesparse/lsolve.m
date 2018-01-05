function x = lsolve (L,b)						    %#ok
%LDLSOLVE solve LL'x=b using a sparse LL' factorization
%
%   Example:
%   x = lsolve (L,b)
%
%   solves the system L*L'*x=b for x.  This is equivalent to
%
%   x = L' \  \ (L \ b) ;
%
%   L is from lchol
%
%   See also LCHOL

%   Copyright 2006-2017, Timothy A. Davis, http://www.suitesparse.com

error ('lsolve mexFunction not found') ;
