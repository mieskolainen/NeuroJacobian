% 2D-camel function to be integrated
% mikael.mieskolainen@cern.ch, 08/2018
function y = f_func(x)

x = x(:);
y = exp(-norm(x-ones(2,1)*0.3).^2/0.2^2) + exp(-norm(x-ones(2,1)*0.7).^2/0.2^2);


end
