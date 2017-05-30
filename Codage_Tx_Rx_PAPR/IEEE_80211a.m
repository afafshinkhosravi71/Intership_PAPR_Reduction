%========================================================================
%          Programme de simulation de la chaine 802.11a
%========================================================================
clear
close all;

%%----------------- Parametres de la sequence OFDM ---------------------%% 

NFFTSize = 64;                          % debit = 6 Mbits/s;
NbSym_Mod = 16;                         % Choix du nombres d'etats de la modulation
NbBit_Mod = nextpow2(NbSym_Mod);        % NbBit_Mod = 4(16), 2(4), 1(2)
NBit = NbBit_Mod * 2^10;                % NBit du signal transmis
Msg_bin = (rand(1,NBit)' > 0.5)*1 ;     % sequence binaire de 0 et 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Modulation du flux binaire ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sg_Mod = modulation (NbSym_Mod, Msg_bin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Modulation OFDM IEEE 802.11a -------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sg_OFDM, NSymb] = Allocation_OFDM(Sg_Mod, NFFTSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- IFFT du signal --------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Application de la iFFT
Sg_OFDM2 = QiFFT(Sg_OFDM, NSymb, NFFTSize);

% Constelleation
scatterplot(Sg_OFDM2);
grid on;
title('Constellation apres modulation OFDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------ Ajout d'intevalle de garde 16 symboles ---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sg_final2, Nb_GI] = Add_GI(Sg_OFDM2, NFFTSize, NSymb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ Visualisation ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleRate = 20e6;           % Define sample rate of baseband signal (Hz)
figure()
plot((1:length(Sg_final2))/sampleRate, real(Sg_final2))
hold on;
plot((1:length(Sg_final2))/sampleRate, imag(Sg_final2))
hold off;
title('OFDM Symbols in the time domain after GI')
xlabel('Time (s)')
ylabel('Amplitude')
legend('I component','Q component');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----- Calcul de puissance moyenne du signal pour la normalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Power = mean(abs(Sg_final2).^2);
Gain_OFDM = sqrt(Power);        % Gain 
% Normalisation de la puissance du signal (puissance moyenne unitaire)
Sg_final = Sg_final2 / Gain_OFDM;  

% Visualisation
figure()
plot((1:length(Sg_final))/sampleRate, real(Sg_final))
hold on;
plot((1:length(Sg_final))/sampleRate, imag(Sg_final))
hold off;
title('Normalised OFDM Symbols in the time domain')
xlabel('Time (s)')
ylabel('Amplitude')
legend('I component','Q component');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

% visualisation
figure();
fsMHz = 20;
st = resample(Sg_final, 2, 1);
[Pxx,W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
xlabel('frequency (MHz)')
ylabel('Power Spectral Density (dB)')
title('Frequency spectrum OFDM signal (802.11a) after GI');

%%--------- Peak Power and average power Calculation ------------------%%
Avg_Power(1:length(Sg_final)) = mean(abs(Sg_final).^2);  
Peak_Power(1:length(Sg_final)) = max(abs(Sg_final).^2);

figure(); hold on;
plot((1:length(Sg_final))/sampleRate, abs(Sg_final).^2,':r','linewidth',2);
plot((1:length(Sg_final))/sampleRate, Avg_Power',':black','linewidth',2);
plot((1:length(Sg_final))/sampleRate, Peak_Power',':green','linewidth',2);
title('OFDM Symbols in the time domain')
xlabel('Time (s)')
ylabel('Power signal')
legend('Envelope signal','Average Power','Peak Power');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ Clipping method -------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_clip = 1.65;                  %% Seuil d'ecretage
Sg_OFDM_clp = clipping(Sg_final, A_clip);

% Constelleation
scatterplot(Sg_OFDM_clp);
grid on;
title('Constellation apres clipping');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------- Calcul PAPR ------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 80 subcarriers per symbol OFDM (64 + 16 IG)
NSubcarriers = NFFTSize + Nb_GI;   

% Dispose les données en matrice (NFFTSize, NSymb) (64 , 22)
Sg_final2 = reshape(Sg_final, NSubcarriers, NSymb); 

PaprSymb_dB = calculPAPR(Sg_final2, NSymb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------- Simulation avec la methode TR ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter_max = 10;                  % Nbre d'iteration

%%------ Application de la methode du Hessien -----------------%
[Sg_OFDM_TR] = Hessien(Sg_OFDM, A_clip, iter_max, NSymb, Gain_OFDM);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- IFFT du signal --------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Application de la iFFT
Sg_OFDM_TR2 = QiFFT(Sg_OFDM_TR, NSymb, NFFTSize);

% Constelleation
scatterplot(Sg_OFDM_TR2);
grid on;
title('Constellation apres TR');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------ Ajout d'intevalle de garde 16 symboles ---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sg_final_TR] = Add_GI(Sg_OFDM_TR2, NFFTSize, NSymb);
Sg_final_clp = Sg_OFDM_clp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------- Visualisation-------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot((1:length(Sg_final_TR))/sampleRate, real(Sg_final_TR))
hold on;
plot((1:length(Sg_final_TR))/sampleRate, imag(Sg_final_TR))
hold off;
title('OFDM Symbols in the time domain after TR')
xlabel('Time (s)')
ylabel('Amplitude')
legend('I component','Q component');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

figure();
fsMHz = 20;
st = resample(Sg_final, 2, 1);
[Pxx,W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
hold on;
st = resample(Sg_final_TR, 2, 1);
[Pxx, W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
hold on;
st = resample(Sg_final_clp, 2, 1);
[Pxx, W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
legend('Spectre de frequence without TR', 'Spectre de frequence with TR', 'Spectre de frequence with clipping');
xlabel('frequency (MHz)')
ylabel('Power Spectral Density (dB)')
title('Frequency spectrum OFDM signal (802.11a) without/with TR');

% Visualisation without/with clipping
figure()
plot((1:length(Sg_final))/sampleRate, real(Sg_final))
hold on;
plot((1:length(Sg_final_clp))/sampleRate, real(Sg_final_clp))
title('Signal with clipping')
xlabel('Time (s)')
ylabel('Amplitude')
legend('I without Clipping','I with Clipping');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

% Visualisation without/with TR 
figure()
plot((1:length(Sg_final))/sampleRate, real(Sg_final))
hold on;
plot((1:length(Sg_final_TR))/sampleRate, real(Sg_final_TR))
hold off;
title('Signal without/with TR method')
xlabel('Time (s)')
ylabel('Amplitude')
legend('I without TR','I with TR');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

%%--------- Peak Power and average power Calculation ------------------%%
Avg_Power_TR(1:length(Sg_final_TR)) = mean(abs(Sg_final_TR).^2);  
Peak_Power_TR(1:length(Sg_final_TR)) = max(abs(Sg_final_TR).^2);

figure(); hold on;
plot((1:length(Sg_final_TR))/sampleRate, abs(Sg_final_TR).^2,':r','linewidth',2);
plot((1:length(Sg_final))/sampleRate, Avg_Power',':black','linewidth',2);
plot((1:length(Sg_final_TR))/sampleRate, Avg_Power_TR','-blue','linewidth',2);
plot((1:length(Sg_final))/sampleRate, Peak_Power',':green','linewidth',2);
plot((1:length(Sg_final_TR))/sampleRate, Peak_Power_TR','-yellow','linewidth',2);
title('OFDM Symbols in the time domain after TR')
xlabel('Time (s)')
ylabel('Power signal')
legend('Power of signal','Average Power without TR','Average Power with TR', 'Peak Power without TR', 'Peak Power with TR');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------- Calcul PAPR ------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Dispose les données en matrice (NFFTSize, NSymb) (64 , 22)
Sg_final_TR2 = reshape(Sg_final_TR, NSubcarriers, NSymb); 
Sg_final_clp2 = reshape(Sg_final_clp, NSubcarriers, NSymb); 

PaprSymbTR_dB = calculPAPR(Sg_final_TR2, NSymb);
PaprSymbclp_dB = calculPAPR(Sg_final_clp2, NSymb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------- visualisation de la CCDF ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
[n, x] = hist(PaprSymb_dB,[0:0.005:15]); 
semilogy(x, 1-cumsum(n)/NSymb,'LineWidth',2, 'color', 'blue'); 
hold on;
[n1, x1] = hist(PaprSymbTR_dB,[0:0.005:15]); 
semilogy(x1, 1-cumsum(n1)/NSymb,'LineWidth',2, 'color', 'red');
hold on;
[n2, x2] = hist(PaprSymbclp_dB,[0:0.005:15]); 
semilogy(x2, 1-cumsum(n2)/NSymb,'LineWidth',2, 'color', 'green');
grid on;
xlabel('PAPR_0 (dB)')
ylabel('Probability(PAPR > PAPR_0)')
title('CCDF du PAPR avec transmission IEEE 802.11a ')
legend('signal without TR', 'signal with TR','signal with Clipping'); 

%%---- Calcul de l'augmentation de la puissance caus� par la TR en dB
Power_Mean_Square = mean(abs(Sg_final).^2); 
Power_Mean_Square_TR = mean(abs(Sg_final_TR).^2); 
Power_increase_dB = 10 * log10(Power_Mean_Square_TR  / Power_Mean_Square);
Power_increase_dB
