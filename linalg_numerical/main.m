% Asm1_19AE10037

clear all

A = [0 4.4 3 -6.6; 0.4 3.6 0 8.4; -2 -6.2 5 0; 1 0 -7.6 3];
b = [-4.65; 4.62; -4.35; 5.97];

gaussian.jordan_eliminate(A, b)