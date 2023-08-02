close all; clear; clc;

[M, rhs, tsol] = import_vectors("vectors.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("operators.csv");

M = spdiags(M,0,length(M),length(M));
% 
% D_0 = reshape(D_0,[sqrt(length(D_0)),sqrt(length(D_0))])';
% D_1 = reshape(D_1,[sqrt(length(D_1)),sqrt(length(D_1))])';
% 
% L_0 = reshape(L_0,[sqrt(length(L_0)),sqrt(length(L_0))])';
% L_1 = reshape(L_1,[sqrt(length(L_1)),sqrt(length(L_1))])';

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';

eig(G_0)

eig(G_1)

eig(G_0'*M*G_0 + G_1'*M*G_1)