
clear all
close all
clc

% The directory of the file to be sent
currentFolder = pwd;
cd (currentFolder);

% Delete this file if it exists because it might cause some problems 
% if exist('Sent_Frame.it', 'file')==2
%   delete('Sent_Frame.it');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------- Transmitter ----------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excuting the it++ transsmision chain
% command='./Transmitter';
% [status,cmdout] = system(command,'-echo');
% disp('Done -->> IQ data is ready ');

% Loading IQ data
clear all

itload Sent_rame.it;               % loading the IQ file
IQData_c = PPDU_Frame_Final(1:end);

clear PPDU_Frame_Final;

itload Sent_Fame.it;               % loading the IQ file
IQData_TR = PPDU_Frame_Final(1:end);

clear PPDU_Frame_Final;

itload Sent_Frame.it;               % loading the IQ file
IQData = PPDU_Frame_Final(1:end);


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
Peak_Power(1:length(idealSpacedPulse)) = max(abs(idealPulse).^2);
Avg_Power(1:length(idealSpacedPulse)) = mean(abs(idealPulse).^2);
% PAPR_dB = 10 * log10(Peak_Power./Avg_Power);

% PAPR's calculation
nSymbol = length(IQData)/80;    % length(IQData)/(nIFFT + nIG);
[papr_dB] = Calcul_papr(IQData, nSymbol);
[papr_dB_cl] = Calcul_papr(IQData_c, nSymbol);
[papr_dB_TR] = Calcul_papr(IQData_TR, nSymbol);

% plot
if (plot_data==1)

    figure(1); hold on;
    plot((1:length(idealSpacedPulse))/sampleRate, real(idealSpacedPulse),'b');
    plot((1:length(idealSpacedPulse))/sampleRate, imag(idealSpacedPulse),'g');
    title('OFDM Symbols in the time domain')
    xlabel('Time (s)')
    ylabel('Amplitude')
    legend('I component','Q component');
    axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])
    
    figure(2); hold on;
    plot((1:length(idealSpacedPulse))/sampleRate, abs(idealSpacedPulse),':r','linewidth',2);
    plot((1:length(idealSpacedPulse))/sampleRate, Avg_Power',':black','linewidth',2);
    plot((1:length(idealSpacedPulse))/sampleRate, Peak_Power',':yellow','linewidth',2);
    title('OFDM Symbols in the time domain')
    xlabel('Time (s)')
    ylabel('Power signal')
    legend('Envelope signal','Average Power','Peak Power');
    axis tight; axisLimits = axis; axis([axisLimits(1:2) 1.2*(axisLimits(3:4))])
    
    figure(3);
    fsMHz = 20;
    % st = resample(idealSpacedPulse,2, 1);
    [Pxx,W] = pwelch(idealSpacedPulse,[],[],4096,20);    
    plot([-2048:2047]*fsMHz/4096,10*log10(fftshift(Pxx))-max(10*log10(fftshift(Pxx))));
%     freq_mask = [-30 -20 -11 -9 9 11 20 30];
%     mask = [-40 -28 -20 0 0 -20 -28 -40];
%     hold on
%     plot(freq_mask, mask)
%     legend('Spectre de frequence', 'Mask de spectre');
    xlabel('frequency (MHz)')
    ylabel('Power Spectral Density (dB)')
    title('Frequency spectrum OFDM signal (802.11a)');
    
    figure(4);
    [n, x ] = hist(papr_dB,(0:0.005:length(papr_dB))); 
    semilogy(x, 1-cumsum(n)/nSymbol,'LineWidth',2, 'color', 'blue');
    hold on;
    [n1, x1 ] = hist(papr_dB_cl,(0:0.005:length(papr_dB_cl))); 
    semilogy(x1, 1-cumsum(n1)/nSymbol,'LineWidth',2, 'color', 'red'); 
    hold on;
    [n11, x11] = hist(papr_dB_TR,(0:0.005:length(papr_dB_TR))); 
    semilogy(x11, 1-cumsum(n11)/nSymbol,'LineWidth',2, 'color', 'green', 'LineStyle', '--');
    xlabel('PAPR_0 (dB)');
    ylabel('CCDF : Probability (PAPR > PAPR_0)');
    title('CCDF du PAPR avec transmission IEEE 802.11a'); grid on;
    legend('Original Signal ', 'Signal with clipping', 'signal with TR');
    
    figure(5);
    Org = imread('lena.jpg');
    Rec = imread('test.jpg');
    [mse, mae, SNR, PSNR] = evaluate (Org, Rec);
    subplot(121); imshow(Org); title('Sented Image');
    subplot(122); imshow(Rec); title('received Image');
    
end
disp('Done -->> IQ data is loaded to the MXG');

itload('modulated_symbols.it')
scatterplot(modulated_symbols)
title('Constellation')

figure()
stem(real(modulated_symbols))
xlabel('Samples');
ylabel('Amplitude');
title('Modulated signal')

gain = sqrt(mean(abs(modulated_symbols_pilots).^2));

itload('Pilot_insertion.it')
figure(); hold on
plot(abs(modulated_symbols_pilots/gain))
xlabel('Samples');
ylabel('Amplitude');
title('Pilot insertion')

itload('data_ofdm_transmited.it')
figure()
plot(real(data_ofdm))
xlabel('Samples');
ylabel('Amplitude');
title('Data ofdm transmited')
scatterplot(fft(data_ofdm))
itload('Clipping.it')

itload('Tone_reservation.it')
t = length(data_ofdm);
figure(); hold on;
plot(abs(data_ofdm(1:t))); 
plot(abs(Tone_reserv(1:length(Tone_reserv))))
plot(A_clip*ones(t), 'g');
xlabel('Samples');
ylabel('Amplitude');
title('TR applied on OFDM signal')
legend('Signal without TR', 'Signal with TR', 'Threshold Aclip = 1.65');


t = length(data_ofdm);
figure(); hold on;
plot(abs(data_ofdm(1:t))); 
plot(abs(Clipping(1:t)))
plot(A_clip*ones(t), 'g');
xlabel('Samples');
ylabel('Amplitude');
title('Clipping applied on OFDM signal')
legend('Signal without clipping', 'Signal with clipping', 'Threshold Aclip = 1.65');

scatterplot(Clipping)
title('constellation afetr clipping');

figure();
fs = 20;
st = resample(data_ofdm,2, 1);
[Px,W] = pwelch(st,[],[],4096,20);    
plot([-2048:2047]*fs/4096,10*log10(fftshift(Px))-max(10*log10(fftshift(Px))));
xlabel('Frequence (MHz)');
ylabel('DSP (dB)');
title('Frequency spectrum OFDM signal');

itload('demod.it')
scatterplot(ofdm_demodulated_symb)
title('constellation à la reception');

