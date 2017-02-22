function L2 = Chol_insert(L, v, c)
% function L2 = Chol_insert(L, v, c)
% Chollesky decomposition update when we insert one column and one row
% suppose that A = L*L' (A is symmetric PSD), L is low triangular
% and A2 = [A v;
%			v'c]
% suppose that A2 = L2*L2', calculate L2 efficienetly 
% Ref: http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf
% Lw = v 
% L2 = [L 0;
%		w' sqrt(c - w'w)]
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 8:26:23 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

%% test mode 
if nargin == 0 
	B = rand(10, 5); 
	b = rand(10, 1);
	B2 = [B, b];
	A = B'*B ;
	A2 = B2'*B2;
	L = chol(A)';
	L3 = chol(A2)';
	%
	v = A2(1:end-1, end);
	c = A2(end, end);
end 
%% main code starts here
w = forsub(L, v);
c2 = sqrt(c - w'*w);
L2 = [L zeros(size(L,1), 1);...
		w' c2];
%% test mode
if nargin == 0 
	disp(norm(L2 - L3));
    L2 = [];
end 


