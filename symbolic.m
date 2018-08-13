% Symbolic model transfer functions and Jacobian/Inverse Jacobian determinants
% mikael.mieskolainen@cern.ch, 08/2018
clear; close all;

addpath ./src

N = 2;  % Input/output dimension
M = 3;  % Hidden layer dimension

z  = sym('z_',  [N 1]);
G0 = sym('G0_', [M N]);
G1 = sym('G1_', [N M]);
b0 = sym('b0_', [M 1]);
b1 = sym('b1_', [N 1]);


%% Code generation:: create un-vectorizing function

filename = sprintf('./src/vec2matN_%d.m', N);
file = fopen(filename, 'w');
k = 1;

fprintf(file, '%% Un-vectorizing function, mikael.mieskolainen@cern.ch \n');
fprintf(file, 'function [G0,G1,b0,b1,dim,biasindex] = vec2matN_%d(x) \n', N);

k = printsymbolic(G0, 'G0', k, file);
k = printsymbolic(G1, 'G1', k, file);
k = printsymbolic(b0, 'b0', k, file);
k = printsymbolic(b1, 'b1', k, file);

fprintf(file, 'dim = %d; \n', k-1);
biasindex = k-(length(b0) + length(b1));
fprintf(file, 'biasindex = %d;\n', biasindex);
fprintf(file, 'end \n');
fclose(file);


%% Analytic transfer functions

close all;

vecsigmoid = @(x) ones(length(x),1) ./ (1 + exp(-x));
sigmoid  = @(x) 1./(1+exp(-x));
modtanh  = @(x) (3/4)*tanh(x) + (1/4)*x;
elu      = @(x) (exp(x) - 1).*heaviside(-x) + heaviside(x) .*x;
selu     = @(x) (exp(x) - 1).*heaviside(-x) - heaviside(x) .* (exp(-x) - 1);
melu     = @(x) (exp(x) - 1).*heaviside(-x) + heaviside(x) .* (exp(-x) - 1) + 1;
relu     = @(x) heaviside(x).*x;
softplus = @(x) log(1 + exp(x)); 

y = G0*z + b0;     % Input transfer layer
y = tanh(y);        % Hidden response function
y = G1*y + b1;     % Transfer layer
g = vecsigmoid(y); % Output response function

% Create model function
matlabFunction(g, 'File', sprintf('./src/gN_%d', N), 'Vars', {z, G0, G1, b0, b1});


%% Plot transfer functions
close all;

x = linspace(-3,3,1e3);

f = figure;
plot(x, sigmoid(x)); hold on;
plot(x, tanh(x), 'linewidth', 1.2);
plot(x, elu(x));
plot(x, selu(x), '-.');
plot(x, melu(x), '-');
plot(x, relu(x), '--');
plot(x, softplus(x));
plot(x, heaviside(x), 'k:');

axis square;
l = legend('sigmoid','tanh','elu','symmetric elu', ...
    'mirror elu','relu', 'softplus', 'heaviside', ...
    'location','northwest');
set(l,'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');

filename = sprintf('transferfunc');
print(f, sprintf('./figs/%s.pdf', filename), '-dpdf');
cmd = sprintf('pdfcrop --margins 10 ./figs/%s.pdf ./figs/%s.pdf', filename, filename); system(cmd);


%% Jacobian w.r.t to input (z)
% Generating Jacobians can take several hours!

J = jacobian(g, z);

% Generate Jacobian determinant
dJ = det(J);

% Generate inverse Jacobian determinant
dinvJ = det(inv(J));

% Generate function files
matlabFunction(dJ,   'File', sprintf('./src/detJ_N%d', N),    'Vars', {z, G0, G1, b0, b1});
matlabFunction(dinvJ,'File', sprintf('./src/detinvJ_N%d', N), 'Vars', {z, G0, G1, b0, b1});

