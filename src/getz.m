% Input prior distribution p(z) for n-dimensions
% mikael.mieskolainen@cern.ch, 08/2018
function [z,p] = getz(n)

%{
% Uniform [a,b] x ... x [a,b]
a = 0;
b = 1;
z = a + (b-a)*rand(n,1);

% Total pdf value 
p = 1;
for i = 1:n
   p = p * (1/(b-a)); 
end
%}

mu    = 0;
sigma = 1;

% Unit normal (Gaussian) distributed for n-variables (independent)
z = mu + randn(n,1)*sigma;

% pdf value
p = 1;
for i = 1:n
   p = p*normpdf(z(i),mu,sigma); 
end
%}
end