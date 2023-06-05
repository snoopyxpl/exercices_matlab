%G3E_Champagne_Mathis.m 
%Fichier matlab contenant les réponses aux questions A à C de la section 3.7 du livrable
%APP signal (2022)
%Mathis Champagne

clear variables;
clc;
close all;

%% A
% Définition des paramètres
Te = 1/8000; % Période d'échantillonnage
t = 0:Te:10e-3; % Vecteur temps de 0 à 5 ms avec un pas de Te
y = ones(size(t)); % Signal à 1V pour tout t
y(t>5e-3) = 0; % Signal à 0V après 5 ms

% Affichage du signal
figure;
subplot(2,1,1);
plot(t*1000,y);
xlabel('Temps (ms)');
ylabel('Amplitude (V)');
title('Signal y(n)');

% Calcul et affichage de la FFT
Y = fft(y);
%Avec une telle définition de la fonction y, on obtient sur une dimension
%temporelle (donc positive), une demie porte d'amplitude 1V. On doit s'attendre alors
%d'après les calculs théoriques de l'excercice 3.3 à avoir un résultat qui s'approche d'une "moitié" de
%sinus cardinal comme résultat de la fft.
f = (0:length(Y)-1)/(length(Y)*Te);
subplot(2,1,2);
plot(f,abs(Y));
xlim([0 2e3]);
xlabel('Fréquence (Hz)');
ylabel('Amplitude (V)');
title('FFT de y(n)');

%% B

% Définition des paramètres
Fe = 8000; % Fréquence d'échantillonnage
D = 1; % Durée des signaux en secondes
t = 0:1/Fe:D; % Vecteur temps de 0 à D avec un pas de 1/Fe
f1 = 1440; % Fréquence du premier signal
f2 = 2000; % Fréquence du deuxième signal

% Création des signaux sinusoïdaux
x1 = sin(2*pi*f1*t);
x2 = sin(2*pi*f2*t);

% Multiplication des signaux
y = x1 .* x2;

% Affichage des signaux
figure;
subplot(2,1,1);
plot(t,y);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal y = x1 * x2');

% Calcul et affichage de la FFT
Y = fft(y);
f = (0:length(Y)-1)/(length(Y)/Fe);
%𝑥(𝑡) = sin(2𝜋𝑓1𝑡) × sin(2𝜋𝑓2𝑡) = (1/2)cos[2𝜋(𝑓2−𝑓1)𝑡] − (1/2)cos[2𝜋(𝑓2+𝑓1)𝑡]
%On s'attend donc à avoir un pic de fréquence 560 hz(f2-f1) et un pic à
%3440hz. On s'attend également à avoir des pics de fréquences symétriques
%par rapport à Fe/2 qui correspondent aux harmoniques, car on fait une FFT
%sur une multiplication de sinus.
subplot(2,1,2);
plot(f,abs(Y));
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('FFT de y = x1 * x2');

%% C
file = "alarmes.mp3";
[x,Fs] = audioread(file); 
%NoiseSignalBinaryMatrix = detection(file); %appel du sous programme 
%NoiseSignal = x.*NoiseSignalBinaryMatrix;

%Sensibilité du micro (dBV)
S=-48;
%Gain du micro (dB)
G=30;
%P Référence
Pref = -94;
%Puissance SPL (dB SPL)
Seuil_P_spl = 80;
%durée en seconde (seuil supérieur)
Dt =1;
%durée en seconde (seuil inférieur)
dt =0.5;
Indicedt = Fs/2;
%Seuil dbW
Seuil_dbW = Seuil_P_spl + G + S +Pref;
Seuil_dBm = Seuil_dbW +30;
disp(Seuil_dbW);

%calcul de la puissance du signal:
nbr = length(x);
Te = 1/Fs;
t = (0:nbr-1)*Te;
% Durée de la demi fenêtre d'estimation de la puissance
dFen = 0.05; % demi fenetre de 50 millisecondes
% Convertir dFen en nombre d'échantillons K
K = round(dFen * Fs);
% Initialiser la matrice de puissance
len = length(x);
powerMatrixW = zeros(len, 1);
powerMatrixDbm = zeros(len, 1); %puissance en dBw (dBm)
for n = K+1:nbr-K
    %formule de la puissance (en W)
    powerMatrixW(n) = (1 / (2*K+1)) * sum(x(n-K:n+K).^2);
end
% Convertir la puissance en dBm
powerMatrixDbm = 10*log10(powerMatrixW) + 30;

BinaryMatrix = zeros(len,1);
BinaryMatrix(powerMatrixDbm > Seuil_dBm)=1;
Binarychange = diff(BinaryMatrix);
Start = find(Binarychange == -1); %Quand on passe de 0 à 1, début de la plage
End = find(Binarychange == 1); %Quand on passe de 1 à 0, fin de la plage
%disp(Start);

FalseSilenceIndex = Start(End - Start < Indicedt); %indice du début des faux silences
frames = [];
for i = 1:length(FalseSilenceIndex)
    start = FalseSilenceIndex(i);
    stop = End(find(End > start,1));
    range = start:stop;
    if size(range)<Indicedt
        frames = horzcat(frames,range);
    end
end

falseSilenceFrame = zeros(size(BinaryMatrix));
falseSilenceFrame(frames) = 1;

%Bruits
FalseNoiseBinary = or(BinaryMatrix,falseSilenceFrame);
FalseNoiseSignal = x.*FalseNoiseBinary;

Binarychange = diff(FalseNoiseBinary);
Start = find(Binarychange == 1);
End = find(Binarychange == -1);

RealNoiseBinary = zeros(size(FalseNoiseBinary));
for i = 1:length(Start)
    if (End(i) - Start(i)) >= Fs         % Durée du bruit doit etre > 1 seconde
        RealNoiseBinary(Start(i):End(i)) = 1;  % Si detecter comme vrai son qui dure 1s on met 1 binaire
    end
end
%disp(RealNoiseBinary);
RealNoiseSignal = x.*RealNoiseBinary;

figure(1);
plot(t,RealNoiseSignal);
xlabel('Temps (s)');
ylabel('Valeur');
title('Vrais bruits détectés');

Binarychange = diff(RealNoiseBinary);
Start = find(Binarychange == 1);
End = find(Binarychange == -1);

NoiseRange = {};
for i = 1:length(Start)
    start = Start(i);
    stop = End(find(End>start,1));
    range = start:stop;
    NoiseRange{i} = x(range);
    Noise = NoiseRange{i};
    %1.
    % Calcul de la densité spectrale de puissance
    N = length(Noise);
    f = (-N/2:N/2-1)*Fs/N;
    X = fftshift(fft(Noise))/N; %centré en 0 pour une meilleure visibilité
    Pxx = 2*abs(X).^2/Fs;
    % Affichage de la densité spectrale de puissance
    figure;
    plot(f, Pxx);
    xlabel('Fréquence (Hz)');
    ylabel('Densité spectrale de puissance');
    title('Densité spectrale de puissance du signal sonore');
    
    %2.
    %calcul de la puissance du signal à partir des échantillons sonores
    P_time = 1/length(Noise) * sum(Noise.^2);
    % Calcul de la FFT du signal
    X = fft(Noise);
    % Calcul de la puissance à partir des coefficients de Fourier
    P_freq = sum(abs(X).^2);
    % Afficher les deux puissances
    fprintf('Puissance à partir des échantillons sonores : %.2f\n', P_time);
    fprintf('Puissance à partir des coefficients de Fourier : %.2f\n', P_freq);
    
    %3.
    % calcul de la DSP
    N = floor(length(Noise)/2)*2; % nombre d'échantillons (on s'assure d'avoir un nombre paire pour la division, sinon on a une erreur car il nous faut un entier)
    Xdft = fft(Noise); % transformée de Fourier
    Xdft = Xdft(1:N/2+1); % on ne considère que la première moitié des fréquences
    DSP = (1/(Fs*N)) * abs(Xdft).^2; % calcul de la DSP
    DSP(2:end-1) = 2*DSP(2:end-1); % prise en compte de la symétrie de la DSP
    % détection de la fréquence maximale
    [~, idx] = max(DSP); % recherche de l'indice de la valeur maximale
    f_max = (idx-1)*Fs/N; % calcul de la fréquence correspondante
    % affichage du résultat
    disp(['La fréquence correspondant à la DSP maximale est de ', num2str(f_max), ' Hz']);

    %4.
    N = length(Noise); % Longueur du signal
    DSP = abs(fft(Noise)).^2 / N; % Calcul de la DSP
    DSPdb = 10*log10(DSP); % Conversion en dB
    % Identifier la fréquence à -30dB
    f = linspace(0, 1, N)*Fs; % Vecteur de fréquences
    f_30db = find(DSPdb <= -30, 1, 'first'); % Identifier la fréquence à -30dB
    % Déterminer les limites de la bande à -30dB
    f_min = f_30db; % Limite inférieure de la bande
    f_max = find(DSPdb(f_min+1:end) >= -30, 1, 'last') + f_min; % Limite supérieure de la bande
    % Calculer l'énergie contenue dans la bande
    E_band = trapz(f(f_min:f_max), DSP(f_min:f_max)); % Intégrale de la DSP sur la bande à -30dB
    E_total = trapz(f, DSP); % Intégrale de la DSP sur l'ensemble des fréquences
    Pct_E_band = E_band / E_total * 100; % Pourcentage d'énergie contenue dans la bande à -30dB
    % Afficher les résultats
    disp(['La bande spectrale à -30dB s etend de ', num2str(f(f_min)), ' Hz a ', num2str(f(f_max)), ' Hz.']);
    disp(['Le pourcentage d''énergie contenue dans cette bande est de ', num2str(Pct_E_band), '%.']);

    %5.
    %window_size = round(0.25 * Fs);
    %noverlap = 0;
    %nfft = 2^nextpow2(window_size);
    %subplot(2, 1, 2);
    %spectrogram(Noise, window_size, noverlap, nfft, Fs, 'yaxis');
    %title('Spectrogramme');
    %xlabel('Temps (s)');
    %ylabel('Fréquence (Hz)');

    disp('-----------------------------------------------------------------');
end

