% Asm4_19AE10037
% object oriented approach is taken, please refer to gaussian.m for the detailed methods and modules

clear all

A = [0 4.4 3 -6.6; 0.4 3.6 0 8.4; -2 -6.2 5 0; 1 0 -7.6 3];
b = [-4.65; 4.62; -4.35; 5.97];

% Solution to Ax=b from Gaussian Elimination
% x_gaussian = gaussian.gauss_eliminate(A,b) 

% Inverse using Gauss-Jordan Elimination
% A_inv = gaussian.jordan_inverse(A)

% Solution to Ax=b from Gauss-Jordan Elimination
% x_jordan = gaussian.jordan_eliminate(A,b)

% Solution to Ax=b from Jacobi Iterations
x_jacobi = jacobi.simultaneous(A,b)

% x_T = T(A,b,[0 0 0 0]', 1e-5)