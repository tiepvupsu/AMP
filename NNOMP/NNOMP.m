function x = NNOMP(y, Phi, K)
% Non-Negative Orthogonal Matching Pursuit
% from "On the Uniqueness of Nonnegative Sparse 
% Solutions to Underdetermined Systems of Equations"
% by Bruckstein et al.
% y : data
% Phi : sensing matrix
% K : sparsity
% 
% Gunnar Atli Sigurdsson, 2013

[N2,N] = size(Phi);
x = zeros(1,N);
S = []; % positions indexes of components of s
res = y; % first residual
PhiS = []; % Matrix of the columns used to represent y
normPsquared = shiftdim(sum(Phi.^2,1));

for t=1:K;
    e = norm(res)^2 - max(Phi'*res, zeros(N,1)).^2./normPsquared;
    [~,j]=min(e);
    S = [S j];
    PhiS = [PhiS Phi(:,j)];
    % x_est = lsqnonneg(PhiS,y);
    x_est = myNNQP(PhiS, y, []);
            % rS = y - A(:,S)*xS;
    res = y- PhiS*x_est;
    x(S) = x_est;
end;