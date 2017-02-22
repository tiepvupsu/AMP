function x= forsub(L, b)
% function x= forsub(A, b)
% forward substitution. 
% L is a lower triangular matrix
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 8:35:29 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
if nargin == 0 
	d = 30;
	L = tril(rand(d, d));
	b = rand(d, 1);
	tic; 
	x1 = L\b;
	t1 = toc; 
	tic
	x2 = inv(L)*b;
	t2 = toc;
	tic;
end 
opts.LT = true;
x = linsolve(L, b, opts);
if nargin == 0 
	t3 = toc;
    size(x1)
    size(x2)
    size(x)
% 	disp([x1 x2 x]);
	disp([t1 t2 t3]);
    x = [];
    
end 
