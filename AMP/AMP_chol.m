function [x] = AMP_chol(y, A, lambda, Rho, opts)
% function [x] = AMP_chol(y, A, lambda, Rho, opts)
% using cholesky decomposition
% Abar = [A, sqrt(lambda)*I]
% ybar = [y, zeros]
    if nargin == 0 
        clc;
        d = 10;
        k = 50;
        y = normc(rand(d, 1));
        A = normc(rand(d, k));
        Rho = .2*rand(k, 1);% - .01;
        lambda = 0.1;
        opts.pos = 1;
    end 
    %
    if nargin == 4
        opts.pos = 0;
    end
    %%
    function B = mat(A, S)
        B = A(:, S);
    end
    %
    Aty = A'*y;
    function [xS, rS] = update_xr(L, S)
        % solve xS = argmin||yS - ASxS\|_2 where AS'*AS = L*L'
        % solution AS'*AS*xS = AS'*y =Aty(S,1)
        w = Aty(S, 1);
        v = forsub(L, w);
        xS = backsub(L', v);
        rS = y - A(:,S)*xS;
    end 

    %%
    [d, k] = size(A);
    y = [y; zeros(k, 1)];
    A = [A; sqrt(lambda)*eye(k)];
    %%
	k = size(A,2);
    % selected set 
    S = find(Rho < 0);
    % unselected set:
    U = setdiff(1:size(A,2), S);
    % Initial L and rS, xS 
    if numel(S) > 0
        AS = A(:, S);    
        L = chol(AS*AS')';
        [xS, rS] = update_xr(L, S);
    else
        L = [];
        xS = [];
        rS = y;
    end
    %
    lambda2 = 1/(1 + lambda);
    iter = 0;    
    AtA = A'*A;
	while numel(U) > 0        
        iter = iter + 1;
        Sold = S;
		% --------------- check if insert ------------------------- 
        if ~opts.pos       
            improve_insert = Rho(U) - lambda2*(A(:,U)'*rS).^2;
        else 
            improve_insert = Rho(U) - lambda2*max(0, A(:,U)'*rS).^2;
        end 
        [USbar, i] = min(improve_insert);
        % check if remove 
        VSbar = 1; % a positive number 
        if numel(S) > 0
            improve_remove = (1 + lambda)*xS.^2 + ...
                2*(A(:,S)'*rS).*xS - Rho(S);
            [VSbar, j] = min(improve_remove);
        end 
        if min(USbar, VSbar) > 0
            break
        end
        if USbar < VSbar % update 
%             i_insert = i;
            % v = A(:,S)'*A(:, U(i));
            v = AtA(S, U(i));
            c = 1 + lambda;
            L = Chol_insert(L, v, c);
            % update S and U 
            S = [S, U(i)];
            % U = setdif(U, U(i));
            U(i) = [];
        else % remove 
%             j_remove = j;
            L = Chol_remove(L, j);    
            U = [U, S(j)];
            S(j) = [];
        end 
        %% update 
        if ~opts.pos 
            [xS, rS] = update_xr(L, S);
        else 
%             x0 = zeros(k,1);
%             x0(Sold,:) = xS;
            
            xS = myNNQP(A(:,S), y, []);
            rS = y - A(:,S)*xS;
        end 
    end
    x = zeros(k,1);
    x(S) = xS;
    if nargin == 0
        disp(iter);
        x = [];
    end 
end