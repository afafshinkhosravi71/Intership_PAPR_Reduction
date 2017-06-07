%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------- Siso simple --------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

itload('modulated.it')
scatterplot(modulated_symbols)
grid on;
title('Constellation after modulation')

figure()
stem(real(modulated_symbols))
xlabel('Samples');
ylabel('Amplitude');
title('Modulated signal')

itload('Pilot.it')
figure(); hold on
stem(real(modulated_symbols_pilots(1:64)))
xlabel('Samples');
ylabel('Amplitude');
title('Pilot insertion')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------- Simulation avec la methode TR matlab ----------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modulated_symbols_TR = reshape(modulated_symbols_pilots, 64, 22);

iter_max = 10;                  % Nbre d'iteration
A_clip = 1.65;
NSymb = 22;
Gain_OFDM = 1;

%%------ Application de la methode du Hessien -----------------%
[Sg_after_T] = Hessien(modulated_symbols_TR, A_clip, iter_max, NSymb, Gain_OFDM);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- IFFT du signal --------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NFFTSize = 64;
% Application de la iFFT
Sg_after_T = QiFFT(Sg_after_T, NSymb, NFFTSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------ Ajout d'intevalle de garde 16 symboles ---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sg_after_T] = Add_GI(Sg_after_T, NFFTSize, NSymb);

u = 1;
Sg_after_TR = zeros(1,length(Sg_after_T));
for i = 1:size(Sg_after_T, 1)
    for j = 1:size(Sg_after_T, 2)
        Sg_after_TR(u) = Sg_after_T(i,j);
        u = u + 1;
    end
end

figure(); hold on
stem(real(Sg_after_TR))
xlabel('Samples');
ylabel('Amplitude');
title('Signal after TR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ Modulation OFDM & Visualisation ---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

itload('data_ofdm.it')
Sg_final2 = data_ofdm;
sampleRate = 20e6;

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


itload('Clipping.it')
sampleRate = 20e6;

figure()
plot((1:length(Sg_final2))/sampleRate, abs(Sg_final2))
hold on;
plot((1:length(Clipping))/sampleRate, abs(Clipping))
hold off;
title('Signal clipping')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Without clipping','With clipping');
axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])

figure();
fsMHz = 20;
st = resample(Sg_final2, 2, 1);
[Pxx,W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
hold on;
st = resample(Sg_after_TR, 2, 1);
[Pxx, W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
hold on;
st = resample(Clipping, 2, 1);
[Pxx, W] = pwelch(st, [], [], 4096, 20);    
plot([-2048:2047]*fsMHz/4096, 10 * log10(fftshift(Pxx)) - max(10*log10(fftshift(Pxx))));
legend('Spectre de frequence without TR', 'Spectre de frequence with TR', 'Spectre de frequence with clipping');
xlabel('frequency (MHz)')
ylabel('Power Spectral Density (dB)')
title('Frequency spectrum OFDM signal (802.11a) without/with TR');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------- Calcul PAPR ------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
NSubcarriers = 80;
NSymb = 22;

% Dispose les données en matrice (NFFTSize, NSymb) (64 , 22)
Sg_final3 = reshape(Sg_final2, NSubcarriers, NSymb); 
Sg_after_TR2 = reshape(Sg_after_TR, NSubcarriers, NSymb); 
Clipping2 = reshape(Clipping, NSubcarriers, NSymb); 

PaprSymb_dB = calculPAPR(Sg_final3, NSymb);
PaprSymbTR_dB = calculPAPR(Sg_after_TR2, NSymb);
PaprSymbclp_dB = calculPAPR(Clipping2, NSymb);

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