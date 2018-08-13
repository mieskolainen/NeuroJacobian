% Code generation
% mikael.mieskolainen@cern.ch, 08/2018
function k = printsymbolic(A, A_name, k_in, file)

k = k_in;
fprintf(file, '%s = [\n', A_name);
for i = 1:size(A,1)
    for j = 1:size(A,2)
        fprintf(file, ' x(%d)', k); 
        k = k + 1;
    end
    if (i < size(A,1))
        fprintf(file, ';\n');
    else
        fprintf(file, '\n');
    end
end
fprintf(file, ']; \n\n');

end