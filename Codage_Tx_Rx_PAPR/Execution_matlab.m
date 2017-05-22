%////////////////////////////////////////////////////////////////////////%
%//                IEEE 802.11a OFDM PHYSICAL LAYER                    //%
%//                Author : Taha ALWAJEEH                              //%
%//                Date   : JUNE 2014                                  //%
%////////////////////////////////////////////////////////////////////////%

clear all
close all
clc

% The directory of the file to be sent
currentFolder = pwd;
cd (currentFolder);

% Delete this file if it exists because it might cause some problems 
if exist('Sent_Frame.it', 'file')==2
  delete('Sent_Frame.it');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Transmitter ----------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excuting the it++ transsmision chain
command='./Transmitter';
[status,cmdout] = system(command,'-echo');
disp('Done -->> IQ data is ready ');

% Loading IQ data
clear all
itload Sent_Frame.it;               % loading the IQ file
IQData = PPDU_Frame_Final(1:end);

% Time domain representation
figure(1)
plot(real(IQData(1:end)),'b'); hold on;
plot(imag(IQData(1:end)),'g');
title('Time domain OFDM Symbols')
xlabel('Time (s)')
ylabel('Magnitude')
legend('I (Real) component','Q (Imag) component');

% Time domain representation
% figure(1)
% subplot(211)
% plot(real(IQData(1:7442)),'b'); hold on;
% plot(imag(IQData(1:7442)),'g');
% subplot(212)
% plot(real(IQData(7443:2*7442)),'b'); hold on;
% plot(imag(IQData(7443:2*7442)),'g');
% title('Time domain OFDM Symbols')
% xlabel('Time (s)')
% ylabel('Magnitude')
% legend('I (Real) component','Q (Imag) component');

%819120

% Transmission Parameters
sampleRate = 20e6;                       % Define sample rate of baseband signal (Hz)
Nch = 13;                                % Channel Number Nch = 1, 2, 13 (802.11a)
centerFrequency = (2407 + 5*Nch)*1e6;    % Center Frequency
transmissionPowerdB = -20;               % Transmission Power in dBm
plot_data = 1;                           % variable to plot the data being sent 1-> plot 0-> don't plot

% max size of the IQ vector is > 8e6 
% Normalize amplitude.
scale = max(max(abs(real(IQData))), max(abs(imag(IQData))));
idealPulse = IQData / scale;
idealSpacedPulse = idealPulse;

% Peak Power and average power Calculation
Peak_Power(1:length(idealSpacedPulse)) = (max(abs(idealPulse).^2));
Avg_Power(1:length(idealSpacedPulse)) = (mean(abs(idealPulse)));
PAPR_dB = 10 * log10(Peak_Power./Avg_Power);

% PAPR's calculation
nSymbol = length(IQData)/80;    % length(IQData)/(nIFFT + nIG);
[paprSymboldB] = Calcul_papr(IQData, nSymbol);

% plot
if (plot_data==1)
    figure(2);
    fsMHz = 20;
    % st = resample(idealSpacedPulse,2, 1);
    [Pxx,W] = pwelch(idealSpacedPulse,[],[],4096,20);    
    plot([-2048:2047]*fsMHz/4096,10*log10(fftshift(Pxx))-max(10*log10(fftshift(Pxx))));
    freq_mask = [-30 -20 -11 -9 9 11 20 30];
    mask = [-40 -28 -20 0 0 -20 -28 -40];
    hold on
    plot(freq_mask, mask)
    legend('Spectre de frequence', 'Mask de spectre');
    xlabel('frequency (MHz)')
    ylabel('Power Spectral Density (dB)')
    title('Frequency spectrum OFDM signal (802.11a)');

    figure(3); hold on;
    plot((1:length(idealSpacedPulse))/sampleRate, real(idealSpacedPulse),'b');
    plot((1:length(idealSpacedPulse))/sampleRate, imag(idealSpacedPulse),'g');
    plot((1:length(idealSpacedPulse))/sampleRate, abs(idealSpacedPulse),':r','linewidth',2);
    plot((1:length(idealSpacedPulse))/sampleRate, Avg_Power',':black','linewidth',2);
    plot((1:length(idealSpacedPulse))/sampleRate, Peak_Power','yellow','linewidth',2);
    title('OFDM Symbols in the time domain')
    xlabel('Time (s)')
    ylabel('Amplitude')
    legend('I component','Q component','Envelope','Average Power','Peak Power');
    axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])
    
    figure(4)
    [n_sTR, x_sTR ] = hist(paprSymboldB,(0:0.005:length(paprSymboldB))); 
    semilogy(x_sTR, 1-cumsum(n_sTR)/nSymbol,'LineWidth',2, 'color', 'blue'); 
    xlabel('PAPR_0 (dB)');
    ylabel('CCDF : Probability (PAPR > PAPR_0)');
    title('CCDF du PAPR avec transmission IEEE 802.11a'); grid on;
    
    figure(5);
    Org = imread('lena.jpg');
    Rec = imread('test.jpg');
    [mse, mae, SNR, PSNR] = evaluate (Org, Rec);
    subplot(121); imshow(Org); title('Sented Image');
    subplot(122); imshow(Rec); title('received Image');
    
end
disp('Done -->> IQ data is loaded to the MXG');



itload('modulated_symbols.it')
figure()
stem(modulated_symbols(1:64))
itload('Pilot_insertion.it')
figure()
stem(modulated_symbols_pilots(1:64))
itload('data_ofdm_transmited.it')
figure()
stem(data_ofdm(1:64))
itload('Tone_reservation.it')
figure()
stem(Tone_reserv(1:64))
