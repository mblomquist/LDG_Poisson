close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("operators.csv");

M = spdiags(M,0,length(M),length(M)); 

D_0 = reshape(D_0,[sqrt(length(D_0)),sqrt(length(D_0))])';
D_1 = reshape(D_1,[sqrt(length(D_1)),sqrt(length(D_1))])';

L_0 = reshape(L_0,[sqrt(length(L_0)),sqrt(length(L_0))])';
L_1 = reshape(L_1,[sqrt(length(L_1)),sqrt(length(L_1))])';

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';

T_0 = reshape(T_0,[sqrt(length(T_0)),sqrt(length(T_0))])';
T_1 = reshape(T_1,[sqrt(length(T_1)),sqrt(length(T_1))])';

D_0 = sparse(D_0);
D_1 = sparse(D_1);

L_0 = sparse(L_0);
L_1 = sparse(L_1);

G_0 = sparse(G_0);
G_1 = sparse(G_1);

T_0 = sparse(T_0);
T_1 = sparse(T_1);

A = G_0'*M*G_0 + 0*G_1'*M*G_1;

sol = -pinv(full(A))*(rhs);

norm(sol-tsol)

[sol tsol sol-tsol]







