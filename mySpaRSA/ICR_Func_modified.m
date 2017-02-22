function [x_j,Avg_x,Time_ICR,Objective] = ICR_Func_modified(y,A,AtA,varargin)

%%
% ICR version 1.0, January 20, 2015
%
% ICR solves the spike and slab sparse recovery problem
% with applications in compressed sensing, image restoration
% and reconstruction, sparse regression, classification and
% several other problems.
% This function solves the non-convex problem
%
% arg min_x = || y - A x ||_2^2 + lambda * || x ||_2^2 + rho' * gamma
% where gamma is a vector of size x which gamma_i = Indicator(x_i)
%
% using the ICR algorithm, which is described in "ICR: Iterative Convex
% Refinement for Sparse Recovery Using Spike and Slab Priors" by
% H. S. Mousavi, V. Monga, T. D. Tran, IEEE Signal Processing Letters,
% 2015 (to appear).
%
%
% -----------------------------------------------------------------------
% Copyright (2015): Hojjat s. Mousavi, Vishal Monga and Trac D. Tran
%
% ICR is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%
% Please check http://signal.ee.psu.edu/ICR for the latest version
% of the code
%
%  ===== Required inputs =============
%
%  y: 1D vector of observations of size q
%
%  A: A is a q*p matrix (where 1 is the size of y and p is the size of x)
%
%  AtA: is simply A'*A matrix of size p*p
%
%  ===== Optional inputs =============
%
%  'RHO'   = Vector of size p indicating the coeffient of each gamma_i
%            In case Rho is given as a scalar, it will be assumed as a
%            vector of size p with all equal elements. i.e. different
%            coefficients have the same importance in the model.
%            If not given, the default value is 8e-4 * ones(p,1).
%
%  'LAMBDA' = Regularizer parameter for || x ||_2^2 term in the cost.
%             If not give, the defualt value is 2e-4.
%
%  'TOLERANCESTOP' = Stopping threshold as defined in the algorithm.
%                    If not given, the default value is 0.001.
%
%  'MaxiterA' = maximum number of iterations allowed in the main algorithm.
%               If not give, the defualt is 500.
%
%  'INTERNALITER' = In unconstrained ICR it determines the number of
%                   iterations for computing the refined convex function at
%                   each iteration of the main algorithm.
%
%  'INITIALIZATION' = The initial values pecified for mu_0.
%                     If not given, the default value is 0.01*ones(p,1) for
%                     unconstrained ICR and is ones(p,1) for non-negative
%                     ICR.
%
%  'GROUNDTRUTH' = The true (GroundTruth) underlying x is passed in
%                  this argument if available.
%
%  'ALGORITHM' = Unconstrained or not (1 or 0)
%                Specifies thechoice of algorithm: Unconstrained ICR or
%                non-negative ICR.
%                If not given, the default algorithm is unconstrained ICR.
%
%  'VERBOSE'  = work silently (0) or verbosely (1)
%               If not given, the default value is 1 (verbosely)
%
% ===================================================
% ============ Outputs ==============================
%   x = solution of the main algorithm
%
%
%   Objective = sequence of values of the objective function
%
%   ICR_times = CPU time required for the whole algorithm
%
%
%   mses = sequence of MSE values, with respect to True_x,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
% ========================================================


if (nargin-length(varargin)) ~= 3
    error('Wrong number of required parameters');
end

% Set parameters
Atoms = size(A,2);

% Set the defaults for the optional parameters
Rho = 8e-4 * ones(Atoms,1);
Lambda =  2e-4;
TolStop = 0.001;
MaxIter= 500;
Verbose = 1;
InterIter = 30;
Unconstrained = 1;
x0Avalable = 0;

if (rem(length(varargin),2)==1)
    error('Undefined optional paramaeters. Optional parameters should always be in pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RHO'
                Rho = varargin{i+1};
                if isscalar(Rho)
                    Rho = Rho * ones(Atoms,1);
                elseif size(Rho) ~= size(A'*y)
                    error ('Size of input argument RHO does not match the size of each observation')
                end
            case 'LAMBDA'
                Lambda = varargin{i+1};
            case 'TOLERANCESTOP'
                TolStop = varargin{i+1};
            case 'MAXITER'
                MaxIter = varargin{i+1};
            case 'INTERNALITER'
                InterIter = varargin{i+1};
            case 'INITIALIZATION'
                x_initial = varargin{i+1};
            case 'GROUNDTRUTH'
                x0 = varargin{i+1};
                x0Avalable = 1;
            case 'ALGORITHM'
                Unconstrained = varargin{i+1};
                if Unconstrained == 1
                    fprintf(1,'Unconstrained ICR is running:\n');
                elseif Unconstrained == 0
                    fprintf(1,'Non-Negative ICR is running:\n');
                else
                    error('Wrong parameter for choosing the type of algorithm');
                end
            case 'VERBOSE'
                Verbose = varargin{i+1};
            otherwise
                error(['Undefined option: ''' varargin{i} '''']);
        end;
    end;
end

if Unconstrained == 1
    x_initial = 0.01 *ones(Atoms,1);
else
    x_initial = 0.5*ones(Atoms,1);
end

tic



Stopping = 1;
iter = 0;
while Stopping ==1 && iter< MaxIter
    iter = iter+1;
    
    if iter == 1
        x_j = x_initial;
        Avg_x = x_j +eps;
        Objective(iter) = x_j'*(AtA+Lambda*eye(size(AtA)))*x_j- 2*y'*A*x_j+ Rho' * (abs(x_j)>0) +y'*y;
        ObjDiff(1) = 1;
        SoluDiff(1) = 1;
    elseif iter>1
        x_j_old = x_j;
        
        [x_j] = ConvexRefinement(x_j_old,Avg_x,y,A,AtA,Rho,Lambda,InterIter,Unconstrained);
         
        Stopping = (norm(x_j_old - x_j)) > TolStop;
        Avg_x =  ( ( Avg_x* (iter-1) + (x_j) ) / iter  );
        
        gamma =  abs(x_j)./abs(Avg_x);
        MajorElements = find( abs(gamma - 1 ) <0.5 );
        
        Objective(iter) = x_j'*(AtA+Lambda*eye(size(AtA)))*x_j- 2*y'*A*x_j+ Rho' * (abs(x_j)>0) +y'*y;
        ObjDiff(iter) = Objective(iter) - Objective(iter-1);
        SoluDiff(iter) = norm(x_j - x_j_old);
    end
    if x0Avalable ==1
        MSE(iter) = norm(x_j-x0)^2/Atoms;
    end
    if Verbose ==1
        fprintf(1,'n =%3d, obj= %e, Obj difference = %e, Solution difference = %e, SparsityLevel =%3d\n',iter,Objective(iter),ObjDiff(iter),SoluDiff(iter),sum(abs(x_j)>0));
    end
end

% if Unconstrained == 0
%     options = optimset('Algorithm','interior-point-convex','display','off' ); 
%     RedDic = A(:,MajorElements);
%     RedAtA = RedDic' * RedDic ;
%     RedH=1*( RedAtA );
%     Redf= - 2*y'*RedDic;
%     [Redx_jP,Redfval1P] = quadprog(RedH,Redf,[],[],[],[],zeros(1,length(MajorElements)),[],[],options);
%     Objective(iter) = Redfval1P + y'*y +Rho(MajorElements)' * ones(length(MajorElements),1);
%     x_j = zeros(size(A,2),1);
%     x_j(MajorElements) = Redx_jP;
% end

x = zeros(Atoms, 1);
A_major = A(:, MajorElements);

k_major = size(A_major, 2);
ybar = [y; zeros(k_major, 1)];
Abar = [A_major; sqrt(Lambda)*eye(k_major)];
if Unconstrained == 1
    % x_initial = 0.01 *ones(Atoms,1);
    x(MajorElements) = pinv(Abar)*ybar;
else
    % x = 0.5*ones(Atoms,1);
    x(MajorElements) = lsqnonneg(Abar, ybar);
end
x_j = x;

Time_ICR = toc;

if Verbose ==1
    fprintf(1,'\nFinished the ICR algorithm! Results:\n');
    fprintf(1,'Final number of iterations = %d\n',iter);
    if x0Avalable ==1
        fprintf(1,'MSE of final solution = %6.4e\n', MSE(end));
        fprintf(1,'SNR of final solution = %6.4f dB\n', 20*log10(norm(x0,2)/norm(x_j-x0,2)) );
        fprintf(1,'Support Match measure = %4.2f %%\n', 100*sum(~xor(abs(x0)>0,abs(x_j)>0))/Atoms );
    end
    if Unconstrained == 1
        fprintf(1,'Number of non-zero components (Sparsity Level) = %d\n', sum(abs(x_j)>0));
    else
        fprintf(1,'Number of non-zero components (Sparsity Level) = %d\n', length(MajorElements) );
    end
    fprintf(1,'||A x - y ||_2^2 = %2.4e\n', norm(y-A*x_j)^2);
    fprintf(1,'Penalty due to regularizers = %2.4e\n', Lambda*norm(x_j)^2+Rho'*(abs(x_j)>0) );
    fprintf(1,'Objective function = %6.4e\n',Objective(iter));
    
    fprintf(1,'CPU time  = %4.2f seconds\n', Time_ICR);
    fprintf(1,'\n');
end
