function L2 = Chol_remove(L, i)
% function L2 = Chol_remove(L, i)
% Chollesky decomposition update when we remove i-th column and i-th row of A
% suppose that A = L*L' (A is symmetric PSD), L is low triangular
% and A = [A11 a12 A13;
%			a12' a22 a23';
%			A13' a23 A33];
% L =  [L11  0    0;
%		l21 l22   0;
%		L31 l32 L33];
% A2 = [A11 A13;
%		A13' A33];
% suppose that A2 = L2*L2', calculate L2 efficienetly 
% Ref: http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf
% L2 = [L11 L12;
%		L12' L33new]
% where L33new satisfies: 
%	L33newL33new' = L33L33' + l32*L32'
% then L33new can be obtained by using built-in cholupdate function 
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 8/22/2016 8:26:23 AM
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------

%% test mode 
if nargin == 0 
	B = rand(1000, 500); 
    i = 5;
	B2 = [B(:, 1:i-1 ), B(:, i+1:end)];
	A = B'*B ;
	A2 = B2'*B2;
	L = chol(A)';
	L3 = chol(A2)';
    
	%
end 
%% main code starts here
if i == size(L, 1)
	L2 = L(1:end-1, 1:end-1);
else
	L11 = L(1:i-1, 1:i-1);
	L31 = L(i+1:end,1:i-1);
	l32 = L(i+1:end, i);
	L33 = L(i+1:end, i+1:end);
	L33new = cholupdate(L33', l32, '+')';
% 	L33new = chol(L33' + l32*l32')';
	L2 = [L11, zeros(i-1, size(L33new,2));...
		  L31, L33new];
end
%% test mode
if nargin == 0 
	disp(norm(L2 - L3));
%     disp(L2*L2' - A2)
%     disp(L3)
    L2 = [];
end 


