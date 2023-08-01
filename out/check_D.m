close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("operators.csv");

for i = 1:length(rhs)
    if (abs(rhs(i)) < 1.0e-12)
        rhs(i) = 0;
    end
end

for i = 1:length(rhs)
    if (abs(tsol(i)) < 1.0e-12)
        tsol(i) = 0;
    end
end

D_0 = reshape(D_0,[sqrt(length(D_0)),sqrt(length(D_0))])';
D_1 = reshape(D_1,[sqrt(length(D_1)),sqrt(length(D_1))])';


[rhs tsol D_0*tsol]

norm(rhs - D_0*tsol)