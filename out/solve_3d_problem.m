close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");

M = spdiags(M,0,length(M),length(M)); 

[i0,j0,G0] = import_G_1d("G0.csv");
G_0 = sparse(i0+1,j0+1,G0);
clear i0 j0 G0;

[i1,j1,G1] = import_G_1d("G1.csv");
G_1 = sparse(i1+1,j1+1,G1);
clear i1 j1 G1;

[i2,j2,G2] = import_G_1d("G2.csv");
G_2 = sparse(i2+1,j2+1,G2);
clear i2 j2 G2;

A = G_0'*M*G_0 + G_1'*M*G_1 + G_2'*M*G_2;

sol = -pinv(full(A))*(M*rhs);

sqrt((sol-tsol)'*M*(sol-tsol))
