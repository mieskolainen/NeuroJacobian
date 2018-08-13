% Visualize net
% mikael.mieskolainen@cern.ch, 2018
function f = visualize(x,S,f)

global N

if (nargin < 3)
    f = figure;
end
figure(f);

eval(sprintf('[G0,G1,b0,b1] = vec2matN_%d(x);', N));
gval = zeros(N,S);

% Generate samples
for n = 1:S
    
    % Input prior
    [z,p] = getz(N);
    
    % Generative mapping
    eval(sprintf('g = gN_%d(z,G0,G1,b0,b1);', N));
    
    % Save values
    gval(:,n) = g;
end

subplot(1,2,1);
[z,x] = hist3(gval', [50 50]);

imagesc(x{1},x{2},z');
set(gca,'YDir','Normal');
xlabel('$g_1$','interpreter','latex'); ylabel('$g_2$','interpreter','latex');
axis square; colorbar;

subplot(1,2,2);
surf(x{1},x{2},z');
xlabel('$g_1$','interpreter','latex'); ylabel('$g_2$','interpreter','latex');
axis square;
shading interp;
colormap(hot);

end