close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, D_2, L_0, L_1, L_2, G_0, G_1, G_2, T_0, T_1, T_2] = import_operators_3d("operators.csv");

M = spdiags(M,0,length(M),length(M)); 

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';
G_2 = reshape(G_2,[sqrt(length(G_2)),sqrt(length(G_2))])';

G_0 = sparse(G_0);
G_1 = sparse(G_1);
G_2 = sparse(G_2);

A = G_0'*M*G_0 + G_1'*M*G_1 + G_2'*M*G_2;

sol = -pinv(full(A))*(M*rhs);


sqrt((sol-tsol)'*M*(sol-tsol))

norm(sol-tsol)







