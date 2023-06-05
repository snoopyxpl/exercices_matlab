%G3E_Champagne_Mathis.m 
%Fichier matlab contenant les réponses aux questions A à D de la section 5.5 du livrable
%APP signal (2022)
%Mathis Champagne

clear variables;
clc;
close all;
%% A.
file = "/MATLAB Drive/NotesPiano.wav";
%[signal, Fe] = audioread('/MATLAB Drive/notesFlute.wav');
[signal, Fe] = audioread(file);
Te = 1/Fe;					% periode d'échantillonage
LenSig = length(signal);
deltaf = Fe/LenSig;
t = (0:LenSig-1)*Te;		% vecteur temps échantilloné
f = (0:LenSig-1)*deltaf;

fft_sig = abs(fft(signal));

%paramétres du filtre sur filtreDesigner, On a Sol:392Hz, Fa: 369.99Hz, La:
%415.30Hz. On prend delta_f = 2Hz, delta_2f = 22/2 = 12Hz. On a donc:
%fstop1 = 380Hz, Fpas1 = 390Hz, fpas2 = 394Hz, fstop2 = 404Hz.
%On choisit également Astop = 30dB et Apass = 1dB. On otbient un filtre
%d'ordre 1823.

load("/MATLAB Drive/filtre_sol_16khz.mat");
sigfiltre = filter(Num, 1, signal);
fft_filt  = abs(fft(sigfiltre));

PmSignal = max(signal);
PmSigfiltre = max(sigfiltre);
Rapport = PmSigfiltre/PmSignal;
Msge = sprintf("Pm Signal = %f,  Pm Filtre = %f,  rapport = %1.4f", PmSignal, PmSigfiltre, Rapport);
disp(Msge);

figure(1);
subplot(2,1,1);
plot(t, signal);
xlabel('temps (s)');
ylabel('amplitude (V)');
title('Signal Fichier  NotesFlute.wav');

subplot(2,1,2);
plot(f, fft_sig);
xlabel('frequence (Hz)');
ylabel('|X(f)|');

figure(2);
subplot(2,1,1);
plot(t, sigfiltre);
xlabel('temps (s)');
ylabel('amplitude (V)');
title('Signal après filtrage SOL');

subplot(2,1,2);
plot(f, fft_filt);
xlabel('frequence (Hz)');
ylabel('|X(f)|');

%% B.

%Il y a bien un sol présent dans le fichier Accordpiano.wav, On le voit sur
%le plot.

%% C.

%Pour un filtre d'ordre N, on a N multiplications, N-1 decalages, N-1 additions.

%% D.

%Pour baisser l'ordre du filtre sans altérer la qualité, on peut baisser le
%facteur d'atténuation du filtre (Astop sur matlab).
