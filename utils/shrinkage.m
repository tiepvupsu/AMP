
function X = shrinkage(U, opts)
% function X = shrinkage(U, opts)
% Description:
% 	Soft Thresoding function. Solve one of the following problems:
%		U is a matrix of size d x k 
%		opts.lambda is a positive scalar, a column vector of a matrix 
%		1. if opts.lambda is a scalar:
%			 X = argmin_X 0.5*||X - U|| + opts.lambda||X||_1 
%		2. if opts.lambda is a matrix with the same size as U 
%			X = argmin_X 0.5*||X - U|| + ||opts.lambda .* X||_1 
%		3. if opts.lambda is a column vector, suppose W is a matrix whose each 
%			column is opts.lambda, number of columns is the same as number of 
%			columns of U:
%			X = argmin_X 0.5*||X - U|| + ||W.* X||_1 
%		if `opts.pos = true`, then this function solves each of above problems 
%			with positive constraints.
% Inputs: U: double dense matrix d x k 
%		  opts: scalar of struct 
%			+ scalar: then this is opts.lambda in problem 1 
%			+ struct: 
%			  opts.lambda: scalar, column vector of matrix 
%			  opts.pos: positive constraint (default = false)
% Outputs: X: a sparse matrix in d x k
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 6/9/2016 12:00:35 PM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
	if ~isfield(opts, 'lambda')
		lambda = opts;		
	else 
		lambda = opts.lambda;
	end 
	%%
	if ~isfield(opts, 'pos')
		opts.pos = false;
	end 
	%%
	if numel(opts.lambda) > 1 && size(opts.lambda, 2) == 1 % column vector 
		lambda = repmat(lambda, 1, size(U, 2));
	end 
	%%
	if opts.pos 
		X = max(0, U - lambda);
	else 
		X = max(0, U - lambda) + min(0, U + lambda);
	end
	%% sparsify
% 	X = sparse(X);
end
	