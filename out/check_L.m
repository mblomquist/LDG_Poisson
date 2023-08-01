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

L_0 = reshape(L_0,[sqrt(length(L_0)),sqrt(length(L_0))])';
L_1 = reshape(L_1,[sqrt(length(L_1)),sqrt(length(L_1))])';


% [rhs L_0*rhs]

norm(L_0*rhs)
norm(L_1*rhs)