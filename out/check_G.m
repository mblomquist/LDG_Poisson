close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("operators.csv");

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';

gradx = G_0*tsol;
grady = G_1*tsol;

for i = 1:length(gradx)
    if (abs(gradx(i)) < 1.0e-12)
        gradx(i) = 0;
    end
end

for i = 1:length(grady)
    if (abs(grady(i)) < 1.0e-12)
        grady(i) = 0;
    end
end

for i = 1:length(rhs)
    if (abs(rhs(i)) < 1.0e-12)
        rhs(i) = 0;
    end
end

for i = 1:length(tsol)
    if (abs(tsol(i)) < 1.0e-12)
        tsol(i) = 0;
    end
end

% [tsol gradx rhs];

norm(gradx-rhs)