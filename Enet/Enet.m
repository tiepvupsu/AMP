function x = Enet(y, A, lambda2, lambda, opts)
% function x = Enet(y, A, lambda2, lambda)
% Solve the elastic net problem:
% x = arg\min_{x} 1/2 \|y - Ax\|_2^2 + 1/2 lambda2\|x\|_2^2 + lambda\|x\|_1
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 3:22:33 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
if nargin == 0
	clc;
	d = 512;
	k = 256;
	L = 100;
	A = normc(rand(d, k));
	
	Noise_Std = 0.01;
	
	% y = normc(rand(d, 1));
	lambda = Noise_Std;
	opts.pos = 0;
    if ~opts.pos
        x0 = randn(k, 1);
    else
        x0 = rand(k,1);
    end
	x0 = sparsify(x0, L);
	n0 = rand(d, 1)*Noise_Std;
	y = A*x0 + n0;

	lambda2 = 0.01;
	lambda = 0.05;

	Xinit = zeros(size(A, 2), size(y,2));
end 
Xinit = zeros(size(A, 2), size(y,2));
%%
function cost = calc_f(x)
	cost = 0.5*(normF2(y - A*x) + lambda2* x'*x);
end 
function cost = calc_F(X)
	cost = calc_f(X) + lambda*norm1(X);
end 

BB = A'*A + lambda2*eye(size(A,2));
Aty = A'*y;

function g = grad(X)
%     size(BB*X)
%     size(Aty)
	g = BB*X - Aty;
end 

L = max(eig(BB));
% check_grad(@calc_f, @grad, Xinit);
opts.lambda = lambda;
opts.verbose = 0;
opts.max_iter = 200;
x = fista(@grad, Xinit, L, opts, @calc_F);


if nargin == 0 
	mse1 = mse(x, x0);
	sm1 = SM(x, x0);
	sl1 = sum(x ~= 0);
	disp([mse1, sl1, sm1]);
    x = [];
end 
end 
