% Analytic jacobian mapping simulation
%
% First run symbolic.m, then this.
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;

addpath ./src

global N
N = 2;

%% Evaluate the function shape for visuals

x1val = linspace(0,1,1e2);
x2val = linspace(0,1,1e2);

fval  = zeros(length(x1val), length(x2val));
for i = 1:length(x1val)
    for j = 1:length(x2val)
        fval(i,j) = f_func([x1val(i); x2val(j)]);
    end
end

figure;
imagesc(x1val, x2val, fval);
set(gca,'YDir','normal');
axis square; colorbar;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');


%% First learn uniform (flat) output
close all;

global targetmode;
targetmode = 0;    % Set flat target mode

% Starting values for the mapping parameters
eval(sprintf('[~,~,~,~,dim,biasindex] = vec2matN_%d(ones(1000,1));', N));

xhat0 = randn(dim,1);      % Weight terms
xhat0(biasindex:dim) = 0;  % Bias terms (init with zero)

f = figure;
N_iter = 100;
for i = 1:N_iter
    options = optimset('MaxIter', 20);
    [xhat0,cost] = fminsearch(@costfunc, xhat0, options);
    visualize(xhat0, 1e4, f); subplot(1,2,1);
    title(sprintf('Iter %d/%d', i, N_iter), 'interpreter','latex');
    drawnow;
end

%% Stochastic random search + fminsearch

close all; f = figure;

% Fminsearch
bestx = [];
mincost = 1e9;

targetmode = 1;    % Set real target mode
xhat = xhat0;

N_iter = 100;
for i = 1:N_iter
    
    xhat(1:biasindex-1) = xhat(1:biasindex-1) + randn(length(1:biasindex-1),1)*0.1; % some noise
    
    options = optimset('MaxIter', 10);
    [xhat,cost] = fminsearch(@costfunc, xhat, options);
    if (cost < mincost)
        bestx   = xhat;
        mincost = cost;
        
        visualize(bestx, 1e4, f); subplot(1,2,1);
        title(sprintf('Iter %d/%d', i, N_iter), 'interpreter','latex');
        drawnow;
    end
end

% options = optimoptions(@fminunc);
% [xhat,fval,exitflag,output] = fminunc(@costfunc,xhat,options);

xhat = bestx;


%% Optimize the best candidate
 
options = optimset('MaxIter', 50);
[xhat,cost] = fminsearch(@costfunc, xhat, options);


%% Get matrix representation

eval(sprintf('[G0,G1,b0,b1] = vec2matN_%d(xhat);', N));


%% Calculate integral via the generative/Jacobian mapping

sumw  = 0;
sumw2 = 0; 
S = 1e5;

fprintf('\n');

for n = 1:S
    
    % Input prior
    [z,p] = getz(N);
    
    % Generative mapping
    eval(sprintf('g = gN_(z,G0,G1,b0,b1);', N));
    
    % Both are the same, because 1/det(A) = det(A^{-1})
    % eval(sprintf('AbsInvDetJ = 1/abs(detJ_N%d(z,G0,G1,b0,b1));', N));
    eval(sprintf('AbsInvDetJ = abs(detinvJ_N%d(z,G0,G1,b0,b1));', N));
    
    % Importance weight
    w = f_func(g) / (p * AbsInvDetJ);
    
    sumw  = sumw  + w;
    sumw2 = sumw2 + w^2;
end

I     = sumw/S;
I_err = sqrt((sumw2/S - (sumw/S)^2)/S);

fprintf('JACOB:  Integral = %0.5f +- %0.5f (relative = %0.5f) \n', ...
    I, I_err, I_err / I);

% Calculate integral via direct (naive) MC
sumw  = 0;
sumw2 = 0;
for n = 1:S
    w = f_func(rand(N,1));
    sumw  = sumw + w;
    sumw2 = sumw2 + w^2;
end

I = sumw/S;
I_err = sqrt((sumw2/S - (sumw/S)^2)/S);

fprintf('DIRECT: Integral = %0.5f +- %0.5f (relative = %0.5f) \n', ...
    I, I_err, I_err / I);

