function x_res = CPLEX_sol(y, A, lambda, R)
% Use the function cplexmiqp to solve a mixed-integer quadratic programming problem
%
% The MIQP problem solved in this example is
%   Maximize  x1 + 2 x2 + 3 x3 + x4
%             - 0.5 ( 33x1*x1 + 22*x2*x2 + 11*x3*x3 - 12*x1*x2 - 23*x2*x3 )
%   Subject to
%      - x1 +   x2 + x3 + 10 x4 <= 20
%        x1 - 3 x2 + x3         <= 30
%               x2      - 3.5x4  = 0
%   Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3
%        2 <= x4 <= 3
%   Integers
%       x4

% ---------------------------------------------------------------------------
% File: cplexmiqpex.m
% Version 12.4
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2011. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

try
   % Since cplexmiqp solves minimization problems and the problem
   % is a maximization problem, negate the objective
%    d = 32;
%    k = 64;
%    y = normc(rand(d, 1));
%    A = normc(rand(d, k));
%    R = rand(k, 1);
%    lambda = 0.1;
    k = size(A, 2);

   H = [A'*A + lambda*eye(k), zeros(k); zeros(k, 2*k)];
   f = [-A'*y; R];
   
   M = 100;
   Aineq = [eye(k), -M*eye(k); -eye(k), -M*eye(k)];
   bineq = zeros(2*k, 1);
   
   Aeq = [];
   beq = [];
   lb = [];
   ub = [];
   ctype = strcat(repmat('C', 1, k), repmat('B', 1, k));
%    H = [33   6     0    0;
%       6  22    11.5  0;
%       0  11.5  11    0;
%       0   0     0    0];
%    f     = [-1 -2 -3 -1]';
%    
%    Aineq = [-1  1  1 10;
%       1 -3  1  0];
%    bineq = [20  30]';
%    
%    Aeq   = [0  1  0 -3.5];
%    beq   =  0;
%    
%    lb    = [ 0;   0;   0; 0];
%    ub    = [40; inf; inf; 1];
%    ctype = 'CCCBBB';
   
   options = cplexoptimset;
   options.Diagnostics = 'on';
   
   [x, fval, exitflag, output] = cplexmiqp (H, f, Aineq, bineq, Aeq, beq,...
      [], [], [], lb, ub, ctype, [], options);
   
%    fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
%    fprintf ('Solution value = %f \n', fval);
%    disp ('Values =');
%    disp (x');
    x_res = x(1:k);
catch m
   disp(m.message);
end
end