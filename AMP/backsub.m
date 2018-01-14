function x= backsub(U, b)
% function x= backsub(U, b)
% forward substitution. 
% U is a upper triangular matrix
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 8:35:29 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
if nargin == 0 
	d = 3000;
	U = triu(rand(d, d));
	b = rand(d, 1);
	tic; 
	x1 = U\b;
	t1 = toc; 
	tic
	x2 = inv(U)*b;
	t2 = toc;
	tic;
end 
opts.UT = true;
x = linsolve(U, b, opts);
if nargin == 0 
	t3 = toc;
    size(x1)
    size(x2)
    size(x)
% 	disp([x1 x2 x]);
	disp([t1 t2 t3]);
    x = [];
end 
