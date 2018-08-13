% Un-vectorizing function, mikael.mieskolainen@cern.ch 
function [G0,G1,b0,b1,dim,biasindex] = vec2matN_2(x) 
G0 = [
 x(1) x(2);
 x(3) x(4);
 x(5) x(6)
]; 

G1 = [
 x(7) x(8) x(9);
 x(10) x(11) x(12)
]; 

b0 = [
 x(13);
 x(14);
 x(15)
]; 

b1 = [
 x(16);
 x(17)
]; 

dim = 17; 
biasindex = 13;
end 
