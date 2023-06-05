clear variables;
clc;
close all;

%1

f0 = 1e3;
D = 5; % durée en s
Fe = 2e3; % fréquence d'échantillonnage en Hz
A1 = 1;
t = 0 : 1/Fe : D;

x = A1*sin(2*pi*f0*t);

plot(t,x);
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal sinusoïdal x de 1kHz');

T0 = 1/f0;
