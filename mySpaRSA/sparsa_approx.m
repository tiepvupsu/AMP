function x = sparsa_approx(u, Rho, opts)
% solve:
% x = argmin_x (x - u)^2 + Rho.*gamma
% where gamma(i) = 0 if x = 0 
%                = 1 if x ~= 0 
% Rho is now assumingly nonnegative
% opts.pos = pos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0 
    d = 10;
    u = randn(d, 1);
    Rho = rand(d, 1);
    opts.pos = 1;
end 
if nargin < 3
    opts.pos = 0;
end
  
%%    
if min(Rho) < 0
    fprintf('not implemented yet\n');
else 
    if ~opts.pos
        x = (Rho <u.^2).*u;
    else 
        x = max(0, (Rho <u.^2).*u);
    end 
end 
if nargin == 0 
    disp([x u Rho]);
    x = [];
end 