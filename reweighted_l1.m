function x = reweighted_l1(y, A, lambda, opts)
% function x = reweighted_l1(A, lambda)
% solve: 
% argmin \|Wx\|_1 subject to \|y - A*x\|_2^2 <= lambda 
% using reweighted_l1 method http://stanford.edu/~boyd/papers/pdf/rwl1.pdf
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 2:07:22 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
addpath('SPAMS');
if nargin == 0 
	clc;
	d = 256;
	k = 512;
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

end 
[d, k] = size(A);
w = ones(k, 1);
param.mode = 1;
param.lambda = lambda;
param.pos= opts.pos;
param.numThreads = 1;
param.L = size(A, 1);
max_iter = 100;
iter = 0;
eps = 0.0001;
while iter < max_iter
	iter = iter + 1;
	x = mexLassoWeighted(y,A,w,param);
	w = 1./(abs(x) + eps);
end 
x = full(x);
if nargin == 0 
	mse1 = mse(x, x0);
	sm1 = SM(x, x0);
	sl1 = sum(x ~= 0);
	disp([mse1, sl1, sm1]);
    x = [];
end 