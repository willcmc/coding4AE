% Asm6_19AE10037
% object oriented approach is taken, please refer to jacobi.m for the detailed methods and modules

clear all

% A = [0 4.4 3 -6.6; 0.4 3.6 0 8.4; -2 -6.2 5 0; 1 0 -7.6 3];
% b = [-4.65; 4.62; -4.35; 5.97];

A = [5 -2 3 0;-3 9 1 -2;2 -1 -7 1; 4 3 -5 7]
b = [-1; 2; 3; 0.5]

% Solution to Ax=b from Seidel Iterations
x_seidel = jacobi.seidel(A,b)

% Solution to Ax=b from SOR-Seidel Iterations
x_SOR = jacobi.SOR(A, b, 1.4)