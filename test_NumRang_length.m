clc, close all, clear all

n = 2

[A,B,C,D] = butter(n,1,'s');



NL = 200 
axis = [-2 1 -2 2];
h1 = figure('DefaultAxesFontSize',18);
hold on

L = eig(A);
h(1) = plot(real(L),imag(L),'.k','markersize',15);
[h(2),S] = JohnsonAlg(A,NL,axis,'C','-b');
grid on 
box on

