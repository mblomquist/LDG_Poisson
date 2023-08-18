close all; clear; clc;

[M, rhs, tsol] = import_vectors("v_file.csv");
[D_0, D_1, L_0, L_1, G_0, G_1, T_0, T_1] = import_operators("op_file.csv");

M = spdiags(M,0,length(M),length(M)); 

G_0 = reshape(G_0,[sqrt(length(G_0)),sqrt(length(G_0))])';
G_1 = reshape(G_1,[sqrt(length(G_1)),sqrt(length(G_1))])';

G_0 = sparse(G_0);
G_1 = sparse(G_1);

A = G_0'*M*G_0 + G_1'*M*G_1;

sol = pinv(full(A))*(M*rhs);

sqrt((sol-tsol)'*M*(sol-tsol))

I10 = [1, -0.866025, -0.866025, 0.75;
    0, 0.5, 0, -0.433013;
    0, 0, 0.5, -0.433013;
    0, 0, 0, 0.25;
    1, 0.866025, -0.866025, -0.75;
    0, 0.5, 0, -0.433013;
    0, 0, 0.5, 0.433013;
    0, 0, 0, 0.25;
    1, -0.866025, 0.866025, -0.75;
    0, 0.5, 0, 0.433013;
    0, 0, 0.5, -0.433013;
    0, 0, 0, 0.25;
    1, 0.866025, 0.866025, 0.75;
    0, 0.5, 0, 0.433013;
    0, 0, 0.5, 0.433013;
    0, 0, 0, 0.25;];

Gc_0 = .25*I10'*G_0*I10;
Gc_1 = .25*I10'*G_1*I10;
Ac = Gc_0'*Gc_0 + Gc_1'*Gc_1;

r = rhs;
r = I10'*r;
x = pinv(Ac)*r;
x = I10*x
