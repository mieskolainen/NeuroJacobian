% Optimization loop
% mikael.mieskolainen@cern.ch, 08/2018
function cost = costfunc(x)

global targetmode;
global N;

eval(sprintf('[G0,G1,b0,b1] = vec2matN_%d(x);', N));

EPSILON = 1e-9;

% Number of Monte Carlo samples
S = 1e3;

% Generate samples
cost = 0;
for i = 1:S
    
    % Input prior
    [z,p] = getz(N);
    
    % Evaluate the network
    eval(sprintf('g = gN_%d(z,G0,G1,b0,b1);', N));
    
    % log(abs(Jacobian determinant))
    eval(sprintf('dJ = detJ_N%d(z,G0,G1,b0,b1);', N));
    
    % Functional target
    if (targetmode == 0) % Flat output
        fval = 1;
    end
    if (targetmode == 1) % The proper target function
        fval = f_func(g);
    end
    
    % Cost function
    % log Kullback-Leibler divergence
    cost = cost + (log(p) - log(abs(dJ)) - log(max(fval, EPSILON)) );
end
cost = cost / S;


end