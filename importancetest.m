% Importance sampling 1D-test
%
% Test numerically common features with <self-normalized>
% importance sampling, such as estimate variance instability.
%
%
% mikael.mieskolainen@cern.ch, 2018
clear; close all;


%% Classic numerical integral as the reference (~ exact here in 1D)

lowerlimit = 0;
upperlimit = 1;
numint = integral(@func1D,lowerlimit,upperlimit,'AbsTol',1e-8,'RelTol',1e-5);

Neval = round(logspace(3,4,30));
results = zeros(3,size(Neval,3));


%% Importance sampling proposal distribution

mu = 0.8; sigma = 0.3; % importance ansatz f(x)
q  = @(x) 1.0/sqrt(2*pi*sigma^2) * exp(-(x-mu).^2/(2*sigma^2));

% Pre-normalized normalization constant so that \int_a^b dx q(x)=1 
% for pre-normalized
K  = integral(q, 0, 1); 
qn = @(x) q(x) / K;


%% Plot function on sampling scale [0,1]

f2 = figure;
xval = linspace(0,1, 1e3);
plot(xval, func1D(xval)); hold on; plot(xval, q(xval), 'r--');
xlabel('$x$','interpreter','latex');
l = legend('f(x) : function of interest','q(x) : importance function');
set(l,'interpreter','latex');
axis square;


%% Simulation

for k = 1:length(Neval)
    
    N = Neval(k);
    
    %% Naive MC integral
    r = rand(N,1);
    mcint   = sum(func1D(r))/N;
    mcint_e = sqrt( (sum(func1D(r).^2)/N - mcint^2) / N);
    
    
    %% Importance MC sampling integral
    % Basic idea: draw the sample from a proposal
    % distribution q(x). Then re-weight the integral using
    % importance weights so that the correct distribution is obtained.
    %
    % Sample values between [0,1] using q(x)
    % (we know how to sample from q because it is a simple gaussian)
    r = zeros(N,1);
    for i = 1:N
        done = 0;
        while (done == 0)
        r(i) = mu + sigma*randn(1);
            if (r(i) >= 0 && r(i) < 1)
                done = 1;
            end
        end
    end
    
    % Self-normalized importance sampling (N replaced by sum(w))
    w = 1 ./ q(r);
    impmcint_sn   = sum(func1D(r) .* w) / sum(w);
    impmcint_sn_e = sqrt( sum(w.^2 .* (func1D(r) - impmcint_sn).^2) / sum(w)^2);
    
    % Pre-normalized importance sampling (normalization by N)
    w = 1 ./ qn(r);
    impmcint_pn   = sum(func1D(r) .* w) / N;
    impmcint_pn_e = sqrt( (sum( (func1D(r) .* w).^2 )/N - impmcint_pn^2) / N);
    
    
    %%
    fprintf('Integral values with N = %d \n', N);
    fprintf('Numerical           : %0.5f (~ exact) \n', numint);
    fprintf('MC-direct           : %0.5f +- %0.5f (relative %0.3f %%)\n',  mcint, mcint_e, mcint_e / mcint * 100);
    fprintf('MC-self-normalized  : %0.5f +- %0.5f (relative %0.3f %%) \n', impmcint_sn, impmcint_sn_e, impmcint_sn_e / impmcint_sn * 100);
    fprintf('MC-pre-normalized   : %0.5f +- %0.5f (relative %0.3f %%) \n', impmcint_pn, impmcint_pn_e, impmcint_pn_e / impmcint_pn * 100);
    fprintf('\n');
    
    
    %% Save results
    results(:,k) = [mcint; impmcint_sn; impmcint_pn];
    
end


%% Plot sampling evolution trajectory

f1 = figure;
plot(Neval, ones(size(Neval))*mcint, 'k-'); hold on;

colors = {[1 0 0], [0.1 0.5 0], [0 0.5 1]};
for i = 1:3
    plot(Neval, results(i,:)', 'color', colors{i}, 'linewidth', 1.2);
end
set(gca,'xscale','log');
l = legend('True','Direct MC','SN-Importance MC','PN-Importance MC'); set(l,'interpreter','latex');
xlabel('Number of samples $N$','interpreter','latex');
ylabel('Integral estimate $\hat{I}$','interpreter','latex');
axis square;

