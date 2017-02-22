function x_new = sparsify(x, L)
% function x_new = sparsify(x, L)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, Fri 08 Jul 2016 05:30:37 PM EST
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
    x_new = x;
    k = size(x, 1);
    if L < 1 
        L = round(L*k);
    end 

    k2 = k - L;
    x_new(randsample(k, k2)) = 0;
end

