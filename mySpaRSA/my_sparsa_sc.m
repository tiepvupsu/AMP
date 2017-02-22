function [x] = my_sparsa_sc(y, A, lambda, Rho, opts)
% function [x] = my_sparse(y, A, lambda, Rho, opts)
% Solving the ICR cost function using sparse method
    sc = norm(y, 2);
    
    y0 = y;
    y = y/sc;
%     Rho = Rho;
    if nargin == 0 
        clc;
        d = 64;
        k = 32;
        L = 10;
        A = normc(rand(d, k));
        Noise_Std = 0.01;

        % y = normc(rand(d, 1));
        opts.pos = 1;
        if ~opts.pos
            x0 = randn(k, 1);
        else
            x0 = rand(k,1);
        end
        x0 = sparsify(x0, L);
        n0 = rand(d, 1)*Noise_Std;
        y = A*x0 + n0;
        
        Rho = .2*rand(k, 1);% - .01;
        lambda = 0.01;
        
    end 
    %
    if nargin == 4
        opts.pos = 0;
    end
    if ~isfield(opts, 'max_iter')
        opts.max_iter = 300;
    end 
    if ~isfield(opts, 'pos')
        opts.pos = 0;
    end
    %%
    function c = calc_phi(x)
        c = normF2(y - A*x) + lambda*x'*x + Rho'*(x ~= 0);
    end 
    tolP = 1e-6;
    B = 2*(A'*A + lambda*eye(size(A, 2)));
    Aty = 2*A'*y;
    
    function g = grad(x)
        g = B*x - Aty;
    end
    %% choose alpha _t
    function a_t = choose_alpha(x_new, x_old)
        s_t = x_new - x_old;
        r_t = B*(x_new - x_old); % = grad(x_new)- grad(x_old) but faster
                                % since only one matrix-vector mult
        a_t = s_t'*r_t / (r_t'*r_t);
        if a_t < 0
            fprintf('a_t < 0\n');
        end
    end 
    %% for simplicity, choose M = 0 -- monotone 
    M = 0;
    % initial 
    x_old = zeros(size(A, 2), 1);
    eta = 2;
    a_min = 1;
    a_max = 1000;
    iter = 0;
    max_iter = opts.max_iter;
    M = 10;
    Phi = zeros(1, M);
    Phi(end) = calc_phi(x_old);
    max_phi = Phi(end);
    sgma = 0.1;
    while iter < max_iter 
        iter = iter + 1;
        % choose a_t 
        if iter == 1
            a_t = 1;
        end 
        a_t = a_min;
        while 1 > 0 & a_t < a_max
            u_t = x_old - grad(x_old)/a_t;
            x_new = sparsa_approx(u_t, 2*Rho/a_t, opts);
            %% M = 0 
            if calc_phi(x_new) <= max_phi - sgma/2*a_t*normF2(x_new - x_old)
                break;
            end
            a_t = a_t*eta;
%             a_t = 2*min(a_max, 2*choose_alpha(x_new, x_old));
        end 
%         disp(a_t)
        if norm(x_new - x_old)/norm(x_new) < tolP
            break
        end
        Phi = [Phi(2:end), calc_phi(x_new)];
        max_phi = max(Phi);
        x_old = x_new;
    end 
    x = x_new;
%     x = x*sc;    
    AS = A(:, x ~=0);
    x1 = zeros(size(x_new));
    x1(x_new~=0) = (AS'*AS + lambda*eye(size(AS,2)))\(AS'*y0);
    x = x1;
        
    if nargin == 0 
        mse1 = mse(x, x0);
        sm1 = SM(x, x0);
        sl1 = sum(x ~= 0);
        disp([mse1, sl1, sm1]);
%         disp([x x0]);
        x = [];
    end
    
end 
