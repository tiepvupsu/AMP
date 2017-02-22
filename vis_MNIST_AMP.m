
function vis_MNIST_AMP()
    addpath('utils')
    % for n = [1, 2, 4, 6, 7, 8, 9, 10]
    for n = 10
    close all;
    
%     n = 1;
	filename = strcat('data/vis_', num2str(n), '.mat');
	load(filename);
    X_hist = x_hist;
	iters = size(X_hist, 2);
	mse0 = zeros(1, iters);
	cost0 = zeros(1, iters);
    sm0 = zeros(1, iters);
    filename = strcat('res3_', num2str(n), '.gif');
%     figure(1);
    figure('units','normalized','outerposition',[0 0 1 1])
    S_old = 0; 
    ins = 0;
    rev = 0;
    M = iters;
	for i = 1: M
        xi = x_hist(:, i);
        
        %% ground truth
		subplot(2, 3, 1);
		display_network(x0, 1);
        % t1 = sprintf('Ground Truth, SL: %d', numel(find(x0)));
        t1 = sprintf('Ground Truth');
        title(t1);
        
        xlabel(t1);
        
        %% Cost
		subplot(2, 3, 2);
		cost0(i) = cost(y, A, x_hist(:, i), lambda, Rho);
		plot(1:i, cost0(1:i), 'LineWidth',2)
		xlim([0, iters])
		ylim([0, cost0(1)]);
        xlabel('iterations');
        t1 =  sprintf('cost:  %.4f ', cost0(i));
        title(t1);
        
        %% support match
        subplot(2, 3, 6);
		sm0(i) = SM(x_hist(:, i), x0)*i/iters;
		plot(1:i, sm0(1:i), 'LineWidth',2)
		xlim([0, iters])
		ylim([0, 100]);
        t1 = sprintf('Support Match: %.2f %%', sm0(i));
        title(t1);
        xlabel('iterations');
        
        %% recoverd
		subplot(2, 3, 4);        
		display_network(X_hist(:, i), 1);
%         title('recovered');
        t1 = sprintf('Recovered, SL: %d', sl(xi));
        t1 = sprintf('Recovered');
        title(t1);
        xlabel('aa');
        
        
        %% MSE
		subplot(2, 3, 5);
		mse0(i) = mse(x_hist(:, i), x0);
		plot(1:i, mse0(1:i), 'LineWidth',2);
		xlim([0, iters])
		ylim([0, mse0(1)+.01]);
        xlabel('iterations');
        t1 = sprintf('MSE: %.4f', mse0(i));
        title(t1);
		% drawnow();
        
        %% bar 
        S_new = sl(xi);
        if S_new > S_old 
        	ins = ins + 1;
        elseif S_new < S_old
        	rev = rev + 1;
        end 
        S_old = S_new;

        % fprintf('%d = %d - %d\n', sl(xi), ins, rev);

        
        
        
        subplot(2, 3, 3);
        x = [ins rev];
        str = {'+'; '-';};
		h = bar(x, .5);
        title('insertion and removal');
        set(gca,'ytick',[])
        set(h,'FaceColor','b');
        
		ylim([0, iters+15]);
        xlim([0, 3]);
% 

        

		for ii = 1: 2
		    text(ii, x(ii) + 20, [num2str(x(ii))], 'VerticalAlignment', 'top', 'FontSize', 12)
		end 
		set(gca, 'XTickLabel',str, 'XTick',1:numel(str));
%         set(gca,'position',[0 0 1 1],'units','normalized');
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if i == 1;
          imwrite(imind,cm, filename,'gif', 'Loopcount',1);
        elseif i < M
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',.08);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',5);
        end
	end 
    end

end 

function res = cost(y, A, x, lambda, Rho)
	res = cost_spike_slab(y, A, x, lambda, Rho);
end 

function res = sl(x)
	res = numel(find(x));
end 

function res = sm2(x1, x2)
    s1 = find(x1);
    s2 = find(x2);
    res = numel(intersect(s1, s2)) /numel(union(s1, s2));
end 
