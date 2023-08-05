close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("operators.csv");

M = spdiags(M,0,length(M),length(M)); 

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';

G_0 = sparse(G_0);
G_1 = sparse(G_1);

A = G_0'*M*G_0 + G_1'*M*G_1;

sol = -pinv(full(A))*(M*rhs);


sqrt((sol-tsol)'*M*(sol-tsol))

norm(sol-tsol)
