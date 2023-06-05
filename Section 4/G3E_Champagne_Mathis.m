%G3E_Champagne_Mathis.m 
%Fichier matlab contenant les réponses aux questions A à L de la section 4.6 du livrable
%APP signal (2022)
%Mathis Champagne

clear variables;
clc;
close all;

%% A.
% Définition des paramètres du signal
f0 = 440; % Fréquence du sinus en Hz
A = 1; % Amplitude du sinus en V
Te1 = 1/8000; % Période d'échantillonnage en s
duree = 0.1; % Durée du signal en s
t1 = 0:Te1:duree-Te1; % Vecteur temps

% Génération du signal s1
s1 = A*sin(2*pi*f0*t1);

% Calcul de la FFT de s1
M = 2^11; % Nombre de points pour le calcul de la FFT
S1 = fft(s1, M);
module_S1 = abs(S1);

% Affichage du signal s1 dans le domaine temporel et en fréquentiel
subplot(2,1,1);
plot(t1, s1);
title('Signal s1 dans le domaine temporel');
xlabel('Temps (s)');
ylabel('Amplitude (V)');
subplot(2,1,2);
freq = (0:M-1)/(M*Te1); % Vecteur fréquence
plot(freq, module_S1);
title('Module de la FFT de s1');
xlabel('Fréquence (Hz)');
ylabel('Module');

%On retrouve bien un pic de fréquence à 440HZ et un pic de fréquence à Fe -
%440h = 7560 Hz, qui est le symétrique par rapport à Fe/2. Ce deuxième pic est la duplication résultant 
% de l'échantillonage. Ce résultat est bien le résultat attendu lorsque l'on effectue les calculs de la
%transformée de fourier sur un sinus de fréquence f0, échantilloné à la fréquence Fe. (cf exercice 3.4)

%% B.

% Définition des paramètres du signal
f0 = 440; % Fréquence du sinus en Hz
A = 1; % Amplitude du sinus en V
Te2 = 1/500; % Période d'échantillonnage en s
duree = 0.1; % Durée du signal en s
t2 = 0:Te2:duree-Te2; % Vecteur temps

% Génération du signal s2
s2 = A*sin(2*pi*f0*t2);

% Calcul de la FFT de s2
M = 2^11; % Nombre de points pour le calcul de la FFT
S2 = fft(s2, M);
module_S2 = abs(S2);

% Affichage du signal s2 dans le domaine temporel et en fréquentiel
figure;
subplot(2,1,1);
plot(t2, s2);
title('Signal s2 dans le domaine temporel');
xlabel('Temps (s)');
ylabel('Amplitude (V)');
subplot(2,1,2);
freq = (0:M-1)/(M*Te2); % Vecteur fréquence
plot(freq, module_S2);
title('Module de la FFT de s2');
xlabel('Fréquence (Hz)');
ylabel('Module');

% Dans le domaine temporel, on remarque que l'on perd de l'information par rapport au signal
%d'origine. De plus, les pics de fréquences de la fft sont maintenant des
%valeurs absolues d'u sinus cardinal. le spectre de fréquence du signal échantillonné est déformé et 
%certaines fréquences sont confondues avec d'autres. 
% Cela se traduit par la présence de pics de spectre à des fréquences incorrectes, 
% ce qui rend difficile l'interprétation du spectre de fréquence du signal échantillonné.

%% C.
t1 = 0:1/8000:0.1-1/8000; % temps pour s1
t2 = 0:1/500:0.1-1/500; % temps pour s2

s1 = 1*sin(2*pi*440*t1); % signal s1
s2 = 1*sin(2*pi*440*t2); % signal s2

figure;
plot(t1, s1, 'b', t2, s2, 'mo');
xlabel('Temps (s)');
ylabel('Amplitude');
legend('s1', 's2');

%La superposition des deux signaux dans le domaine temporel montre que le signal s2 a une 
%fréquence d'échantillonnage plus basse que s1, ce qui entraîne une résolution temporelle plus faible. 
% Cela se traduit par une représentation plus grossière du signal s2 dans le domaine temporel. 
% On peut voir que le signal s2 suit le même modèle que s1, mais il est moins lisse car il est échantillonné 
% à une fréquence plus basse. On voit bien que l'on perd de l'information par rapport au signal S1. 
% Le critère de shannon n'est pas respecté car Fe_2 = 500 hz <2Fmax_2 = 880 Hz.
%Avec la fréquence f3 pour relier les points de S2, on aurait donc F3 <
%Fe_2, ainsi le critére de Shannon n'est encore pas respecté et l'échantillonage aura une qualité plus moindre 
% que l'échantillonage avec Fe2. 

%% D.
%Si on déséchantillonne le signal S2 par un filtre interpolateur de fréquence de coupure Fe2/2 = 0.25 kHz, 
% On récupère des fréquences parasytes car la fft de S2 n'est pas bonne,
% dû au fait que l'on a pas échantilloné avec le critère de Shannon. Ainsi
% la reconstruction du signal S2 ne sera pas satisfaisant.

%% E.
%Le critère de Shannon est donc une condition nécessaire et suffisante qui
%permet de passer d'un signal continu a un signal discret en évitant les
%phénomènes de repliement. Cela se traduit par le fait que l'on ne perd
%pas d'information sur le signal d'origine. On peut par exemple remarquer
%cette perte d'information a travers l'échantillonage du signal S2, qui est
%le même que le signal S1 mais qui est beaucoup moins satisfaisant.

%% F.
[x, fs] = audioread('Pi_C_96K.wav');

% La fréquence d'échantillonnage est donnée par la variable fs.
%sound(x, fs);
disp(fs);
%fs = 96000Hz

%% G.
% Tracer le spectre du signal audio
L = length(x);
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
freq = (0:L-1)*(fs/L);
% Calcul de l'occupation du spectre
occupancy = sum(abs(Y).^2)/L^2;
disp(['L’occupation du spectre est ', num2str(occupancy)]);

%figure;
plot(freq,abs(Y))
title('Spectre de Pi_C_96K.wav')
xlabel('Fréquence (Hz)')
ylabel('Amplitude')

%% H
% Afficher le signal temporel. 
figure;
k = 6;
x_dec = x(1:k:end);
t_dec = (0:length(x_dec)-1)/fs*k;
plot(t_dec, x_dec);
xlabel('Time (s)');
ylabel('Amplitude');

% Afficher également le spectre. 
N_dec = length(x_dec);
X_dec = fft(x_dec);
f_dec = (0:N_dec-1)*(fs/k/N_dec);

figure;
plot(f_dec, abs(X_dec));
xlabel('Fréquence (Hz)');
ylabel('Amplitude');

% sound(x_dec, fs/k);
%On remarque que le signal est plus aigu et que sa qualité est altérée.
%Dans le domaine fréquentiel, on remarque également des pics de fréquences
%parasite.

%% I.
%En effectuant une décimation, on échantillone le signal audio à une
%fréquence de fs/k. Ainsi, plus la valeur de k est grande, plus la
%fréquence d'échantillonage sera faible, jusqu'à ne plus respecter le
%critère de Shannon.
%La plus grande valeur acceptable pour k dépend donc de la fréquence maximale du signal. 

%% J.
%Avant de sous-échantillonner un signal, il est important de filtrer le signal avec un filtre 
%passe-bas pour éviter le repliement de spectre.

%% k.

[xe, fs] = audioread('BonneJournee.wav');

% Définition du nombre de bits de quantification
b = 8;

% Quantifier le signal sur b bits (dont 1 bit de signe)
xq = round(xe*2^(b-1))/2^(b-1);

% Calculer l'erreur de quantification
eq = xe - xq;

% Calculer le rapport signal à bruit de quantification (SNR)
SNR_b = 10*log10(var(xe)/var(eq));

% Tracer le rapport signal à bruit de quantification en fonction du nombre de bits
SNR_b_values = zeros(1,16);
for b = 1:16
    xq = round(xe*2^(b-1))/2^(b-1);
    eq = xe - xq;
    if var(eq) == 0
        SNR_b_values(b) = NaN;
    else
        SNR_b_values(b) = 10*log10(var(xe)/var(eq));
    end
end
figure;
plot(SNR_b_values);
xlabel('Nombre de bits');
ylabel('SNR (dB)');
title('SNR en fonction du nombre de bits');

% Tracer le rapport signal à bruit de quantification en fonction du débit
SNR_D_values = zeros(1,16);
for b = 1:16
    D = fs*b;
    xq = round(xe*D)/D;
    eq = xe - xq;
    SNR_D_values(b) = 10*log10(var(xe)/var(eq));
end
figure;
plot(SNR_D_values);
xlabel('Débit (bits/s)');
ylabel('SNR (dB)');
title('SNR en fonction du débit');

% Ecouter le signal quantifié
% sound(xq, fs);
%Plus le nombre de bits est important, plus le SNR est grand, c'est à dire
%que la précision augmente et il y a moins d'erreurs d'approximation.

%% L
%Ainsi, à travers cette étude de signaux, on remarque que pour numériser au
%mieux un signal continu, il faut tout d'abord respecter le critère de
%Shannon pour éviter les effets de repliement de spectre et donc éviter au
%maximum la perte d'information sur le sigal d'origine. Ensuite, il faut effectuer un choix 
%judicieux sur le nombre de bits utilisé. Plus le nombre sera élevé, plus
%la précision du signal numérisé sera grande et donc moins il y aura de
%bruits, c'est à dire moins d'erreur sur les échantillons approximés. Il
%faudra cependant, en fonction de notre application, faire attention au
%débit qui peut être assez important. Il faut donc définir un seuil à partir duquel 
%Le rapport signal à bruit parraît acceptable.
