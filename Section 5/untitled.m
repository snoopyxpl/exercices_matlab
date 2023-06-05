%Exercice 5.3
%Si: fréquence 494hz, La: 466hz,Do: 523Hz
%1.Faire un filtre passe bande centré d'une fréquence f0 = 494Hz

%2.Le filtre en réalité laisse passer les fréquence de f0 -delta_f et
%f0+delta_f. Prendre delta_f en fonction de l'écart de fréquence entre les
%notes. (i.e 28Hz et 57Hz).Comme la pente ne peit pas etre verticale, c'est
%une droite avec un coeff directeur, qui commence a partir de delt_f. Cette
%pente delta_2d ne doit pas atteindre les fréquences voisines (La et Do).
%Il faut donc prendre delta_f = 2hz. Il faut donc églameent delta_2f =14hz.
%Fréquence d'écart la + petite /2.
%fstop1 = fo - delta_2f, 
%Fpas1 = f0-delta_f
%Pas2 = fo+delta_f
%Fstop2= f0+delta_2f
%Filtre d'ordre N: N multiplication, N-1 decalage, N-1 addition.
%Amax = 1dB (Apass), Amin: régler en fonction de l'odre du filtre
%(Astop)(entre 30 et 60 dbm?)
%Export le filtre en mat,le load et appeler le filter avec le nom des
%coeffs (de base num)