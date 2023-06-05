
%pdbw = dp SPL + G + SdbV +Pref
function RealNoiseBinary = detection(file)
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
%Seuil dbW
Seuil_dbW = Seuil_P_spl + G + S +Pref;
Seuil_dBm = Seuil_dbW +30;
disp(Seuil_dbW);

curentfile = file;

% Charger les signaux
[x, Fs] = audioread(curentfile);

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


%dt = Fs/2;
FalseSilenceIndex = Start(End - Start < Fs/2); %indice du début des faux silences
frames = [];
for i = 1:length(FalseSilenceIndex)
    start = FalseSilenceIndex(i);
    stop = End(find(End > start,1));
    range = start:stop;
    if size(range)<Fs/2
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
end





