% Test function
% mikael.mieskolainen@cern.ch, 2018

function y = func1D(x)
y = 0.5*normpdf(x, 0.5,0.3) + 0.2*normpdf(x, 0.5,2);
end