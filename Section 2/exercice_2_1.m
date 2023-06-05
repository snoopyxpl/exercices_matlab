clear variables;
clc;
close all;

%1

%paramètres du signal x1
f1 = 1e3;  % fréquence en Hz
A1 = 1;    % amplitude en V
D = 50e-3; % durée en s
Fe = 20e3; % fréquence d'échantillonnage en Hz

%vecteur temps associé à x1 et x2
t = 0 : 1/Fe : D;

%signal sinusoïdal
x1 = A1 * sin(2*pi*f1*t);

% Tracer le signal
plot(t, x1);

%labels
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title('Signal sinusoïdal x1 de 1kHz');

%échelle appropriée pour l'axe des abscisses et l'axe des ordonnées
axis([0 D -A1 A1]);

%2

% Extraction des indices correspondant aux périodes 5 à 10
periode = 1/f1;
periode5 = 5*periode;
periode10 = 10*periode;
ind_periode5 = find(t>=periode5, 1);
ind_periode10 = find(t>=periode10, 1);

% Extraction des données de x1
y = x1(ind_periode5:ind_periode10);

%nouveau vecteur temps pour y
t_y = t(ind_periode5:ind_periode10);

% on trace y sur le même graphe que x1
hold on;
plot(t_y, y);

%3

%paramètres du signal x2
f2 = 1e2;
A2 =0.3;
offset = 2;

%signal sinusoïdal 
x2 = offset + A2 * sin(2*pi*f2*t);

%on trace x2 sur le même graphe
hold on;
plot(t,x2);

%légende
legend('x1','y','x2');

%4

%Le calcul matriciel est la multiplication de deux matrices à l'aide de l'opérateur *. Pour que deux matrices soient compatibles pour la multiplication, le nombre de colonnes de la première doit être égal au nombre de lignes de la deuxième. Le résultat de la multiplication est une nouvelle matrice.
%Le calcul élément par élément est la multiplication terme à terme de deux matrices ou vecteurs à l'aide de l'opérateur .*. Les deux matrices doivent avoir la même taille pour être compatibles. Le résultat est une nouvelle matrice ou vecteur qui a la même taille que les matrices d'origine.

y2 = x1.*x2; 

%En utilisant des opérations matricielles plutôt que des boucles for, on peut souvent accélérer le calcul et réduire le temps d'exécution du code



