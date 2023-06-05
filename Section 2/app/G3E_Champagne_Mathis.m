%G3E_Champagne_Mathis.m 
%Fichier matlab contenant les réponses aux questions A à F de la section 2.6 du livrable
%APP signal (2022)
%Mathis Champagne

clear variables;
clc;
close all;


%% A
% Définition des paramètres du signal
f0 = 1e3;           % fréquence en Hz
A = 2;              % amplitude en V
Fe = 16e3;          % fréquence d'échantillonnage en Hz
Te = 1 / Fe;        % période d'échantillonnage en s
D = 2;              % durée en s

% Calcul du nombre d'échantillons N
N = D * Fe;

% Génération du signal x(nTe)
n = 0:N-1;
x = A * sin(2*pi*f0*n*Te);

% Affichage des 5 premières périodes du signal
T0 = 1 / f0;
n_periods = 5;
n_samples_to_display = n_periods * T0 / Te;

figure;
plot(n(1:n_samples_to_display), x(1:n_samples_to_display));
title('Signal sinusoïdal discret');
xlabel('Temps (s)');
ylabel('Amplitude (V)');



%% B
V_moyenne = mean(x); %valeur moyenne
Px = x.*x;   %puissance instantannée
M_Px = mean(Px); %puissance moyenne
Veff  = sqrt(M_Px); %valeur efficace

disp('Signal x')
Msge = ['X1 : Vmoy= ',num2str(V_moyenne), 'V   Pmoy= ',num2str(M_Px), 'W  Veff= ',num2str(Veff), 'V'];
disp(Msge)



%% C
%1.
x3 = abs(x);
plot(n(1:n_samples_to_display), x3(1:n_samples_to_display))
V_moyenne = mean(x3); %valeur moyenne
Px = x3.*x3;   %puissance instantannée
M_Px = mean(Px); %puissance moyenne
Veff  = sqrt(M_Px); %valeur efficace

disp('Signal x3')
Msge = ['X1 : Vmoy= ',num2str(V_moyenne), 'V   Pmoy= ',num2str(M_Px), 'W  Veff= ',num2str(Veff), 'V'];
disp(Msge)

%2.
%temps d'échantillonnage
Fs = 1000; % Hz
%durée du signal
T = 5; % secondes
% période du signal
T0 = 1; % secondes
%vecteur temps
t = 0:1/Fs:T-1/Fs;
% Initialiser le signal à 0
x4 = zeros(size(t));
% Calcul du signal à partir des conditions de l'énoncé 
for n = 1:length(t)
    if mod(t(n), T0) <= T0/2
        x4(n) = 2;
    else
        x4(n) = 0.5;
    end
end
V_moyenne = mean(x4); %valeur moyenne
Px = x4.*x4;   %puissance instantannée
M_Px = mean(Px); %puissance moyenne
Veff  = sqrt(M_Px); %valeur efficace

disp('Signal x4')
Msge = ['X1 : Vmoy= ',num2str(V_moyenne), 'V   Pmoy= ',num2str(M_Px), 'W  Veff= ',num2str(Veff), 'V'];
disp(Msge)

%% D
NBP = round(T0/Te);
NBP5=NBP*5; % nombre de point sur 5 période du signal
t5p=(1:NBP5)*Te;

B = 1.414;  %B = √(2×10^((𝑃−30)/10)), avec P = 32 DBM, on trouve: B=1.414.

S = pi/3; 
 
y = B * sin((2 * 3.14 * 1 * 10^3 * t5p) + S);
plot(t5p,y);


%% E

Rxx = xcorr(x);
Rxy = xcorr(x, y);
figure;
subplot(2,1,1);
plot(Rxx);
title('Autocorrélation de x');
xlabel('Décalage (n)');
ylabel('Corrélation');
subplot(2,1,2);
imagesc(Rxy);
title('Intercorrélation de x et y');
xlabel('Décalage (n)');
ylabel('Décalage (m)');
colorbar;

%% F

%Variables et calculs

% Charger le signal
[x, Fs] = audioread('MarteauPiqueur01.mp3');
%[x, Fs] = audioread('Jardin01.mp3');
%[x, Fs] = audioread('Jardin02.mp3');
%[x, Fs] = audioread('Ville01.mp3');
nbr = length(x);

Te = 1/Fs;

t = (0:nbr-1)*Te;

% Durée de la demi fenêtre d'estimation de la puissance
dFen = 0.02; % demi fenetre de 20 millisecondes

% Convertion de dFen en nombre d'échantillons K
K = round(dFen * Fs);

% Initialisation de la matrice de puissance
powerMatrix = zeros(length(x), 1);
powerMatrixDbm = zeros(length(x), 1);

% Calculer la puissance du signal
for n = K+1:nbr-K
    % Appliquer la formule de la puissance (en W)
    power = (1 / (2*K+1)) * sum(x(n-K:n+K).^2);

    % Convertir la puissance en dBm
    power_dBm = 10*log10(power) + 30;
    
    % Ajouter la puissance à la matrice de puissance (en W)
    powerMatrix(n) = power;
    powerMatrixDbm(n) = power_dBm;
end

%Affichage des graphiques

figure(1);

subplot(3, 1, 1);
% Plot du signal original
plot(t, x);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal original du marteau piqueur');

subplot(3, 1, 2);
%plot de la puissance (Watt)
plot(t, powerMatrix);
xlabel('Temps (s)');
ylabel('Puissance (Watt)');
title('Puissance du signal');

subplot(3, 1, 3);
%plot de la puissance (DBM)
plot(t, powerMatrixDbm);
xlabel('Temps (s)');
ylabel('Puissance (dBm)');
title('Puissance du signal');
