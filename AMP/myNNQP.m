function x = myNNQP(A, b, x0)
% function x = myNNQP(A, b)
% Solving x = argmin_{x} \|Ax - b\|_2^2 s.t. x >= 0 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 1:01:22 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
if nargin == 0 
	d = 784; k = 1;
	A = rand(d, k);
	b = rand(d, 1);
end 
if nargin <3 | size(x0) == 0
	x0 = zeros(size(A,2), 1);
end
%%
x = x0;
rho = 1;
B = inv(A'*A + rho*eye(size(A,2)));
max_iter = 100;
iter = 0;
u = zeros(size(x));
z = zeros(size(x));
Atb = A'*b;
while iter < max_iter
	iter = iter + 1;
	x = B*(Atb + rho*(z - u));
	z = max(0, x + u);
	u = u + x - z;
end 
x = z;
%%
end 
