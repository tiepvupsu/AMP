function cost = cost_spike_slab(y, A, x, lambda, Rho)
% function cost = cost_spike_slab(y, A, x, lambda, Rho)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, Fri 08 Jul 2016 05:21:50 PM EST
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

    cost = norm(y - A*x,'fro')^2 + lambda*x'*x + Rho'*(x ~= 0);
end 
