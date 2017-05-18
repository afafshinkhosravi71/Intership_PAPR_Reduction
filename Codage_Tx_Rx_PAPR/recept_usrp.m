load test_xt_16qam_12.mat
%x=read_complex_binary('receipt_siso_16qam12_1.raw');
%plot(real(x))
%ginput(2)
%xt=x(ans(1):ans(2));
%plot(real(xt))
%xt=xt(1250:end);

%xt = resampled_IQData(1:int32(length(resampled_IQData)));
sim('Schmidl_Cox.mdl');
%xt = resampled_IQData(1:320);
%sim('paquet_det_canet_et_al.mdl');
[val,ind]=max(yout(1:320));

plot(yout(1:320))

if (ismember(ind,[155:170]))    % if the maximum is around 160 it means that the frame is almost synchronized

    % Frequency offset value 
    %fco=1*yfreq(ind);
    %corr=exp(-j*2*pi*(fco/20e6)*[1:length(resampled_IQData)]);
    %figure(1)
    % plot(real(corr));hold on;
    % plot(imag(corr),'g');hold on;
    % plot(abs(corr),':r');
    %y_CFO = PPDU_Frame_Final.*exp(j*2*pi*0.005*[1:length(PPDU_Frame_Final)]'/64);
    % Frequency offset compensation
    corr=exp(-j*2*pi*(1000/20e6)*[1:length(xt)]); %1106
    xtc=xt.*corr';
    %freq_cor_IQData=resampled_IQData.*corr';
    freq_cor_IQData=xtc;
    % Adjusting the begining of the ofdm frame
    if ((ind-160+1)>0)
        freq_cor_IQData=freq_cor_IQData(ind-160+1:end);
    else
        nbrz=160-ind;
        freq_cor_IQData=[zeros(nbrz,1);freq_cor_IQData];
    end
    
    %plot(real(freq_cor_IQData(1:160)));

else % search for beginning of the frame 
    
   
    
    % search for the peaks which correspond to the begining of the
    % frames the then search for the first peak that correspond to the
    % begininng of our transmitted data
    
    [pks,locs] = findpeaks(yout,'MINPEAKDISTANCE',320);

    % thersholding
    thresh_min = max(pks)/40;
    pks(pks < thresh_min) = 0;

    plot(yout);hold on;
    plot(locs,pks,'k^','markerfacecolor','r')

    indpks = find(pks);
    locspks=locs(indpks);
    
    k=1;
    ind(2)=length(resampled_IQData);
    for(i=1:length(locspks))
        if (i~=length(locspks))
            if (locspks(i+1)-locspks(i)>29500)
                %normally 30000 but a small margin was added
                % this value was added at the transmitter side at the end
                % of the transmitted data to be use for the synch.
                ind(k)=locspks(i+1);
                if (k==2)
                     break;
                end
                k=k+1;
            end
        end
    end
   
    % Frequency offset value 
    fco=1*yfreq(ind(1));
    corr=exp(-j*2*pi*(fco/20e6)*[1:length(resampled_IQData)]);
    figure(3)
    plot(real(corr));hold on;
    plot(imag(corr),'g');hold on;
    plot(abs(corr),':r');
    
    % Frequency offset compensation
    freq_cor_IQData=resampled_IQData.*corr';
    
    % Adjusting the begining of the ofdm frame
    freq_cor_IQData=freq_cor_IQData(ind(1)-160+1:ind(2)-160+1);
    %plot(real(freq_cor_IQData(1:160)));
    if (ind(2)==length(resampled_IQData))
        disp(['Memory limits on the MXA might have caused data acquisition to be truncated, Check !!'])
    end
end

%----------------------------------------------%
% AGC
% This also will be done for every frame in the it++ receiver to get a
% better estimattion

% Maximum value criterion 
% max1=1.1465; 
% This is the max value of the transmitted Short Training Seq.
% max2=max(abs(resampled_IQData(1:160)));
% AGC=max1/max2;
% resampled_IQData=AGC*resampled_IQData;

% Average value criterion
mean1= 0.8663;
% This is the mean value of the transmitted Short Training Seq.
mean2=mean(abs(freq_cor_IQData(1:160)));
AGC=mean1/mean2;
freq_cor_IQData=AGC*freq_cor_IQData;

%----------------------------------------------%
% Save the data to be fed back to the IT++ receiver

% Delete this file if it exists because it might cause some problems 
if exist('IQDATA_matlab.it', 'file')==2
  delete('IQDATA_matlab.it');
end

file_name='IQDATA_matlab.it';       % file name to save the acquired data
disp(['Saving... IQDATA_matlab.it'])
itsave(file_name, freq_cor_IQData);

%plot(real(freq_cor_IQData(1:end)),':r');hold on;
%plot(imag(freq_cor_IQData(1:end)),':b')

%----------------------------------------------%
% Time synchonization and constellation verifications

% Time synchonization verification
synch_verification(freq_cor_IQData);

% constellation verification

% this is just to find how many symbols we should trace for this first 
% channel estimation 

%xt = resampled_IQData(1:int32(length(freq_cor_IQData)));
sim('Schmidl_Cox.mdl');
   
[pks,locs] = findpeaks(yout,'MINPEAKDISTANCE',320);

% thersholding
thresh_min = max(pks)/20;
pks(pks < thresh_min) = 0;
    
plot(yout);hold on;
plot(locs,pks,'k^','markerfacecolor','r')

indpks = find(pks);
locspks=locs(indpks);

nbsym=floor((locspks(2)-locspks(1))/80)-4;% -4 because we dont want to trace the constellation of the first 4 symbols

%nbsymb=28;
plotdata=freq_cor_IQData(locspks(101)-160:locspks(102)-160);
% nbsym=150; %or adjust manually the number of symbols you want to trace
constellation(freq_cor_IQData,nbsym);

% Transfer functions
plot_data=1;          % Frequency Response
transfer_function(freq_cor_IQData,plot_data)
plot_data=2;          % Impulse Response
transfer_function(freq_cor_IQData,plot_data)

%----------------------------------------------%
% IT++ Receiver

% Excuting the it++ reception chain
command='/home/boeglen/Project_Taha_M2_Linux_Windows/System/Receiver_Linux/receiver/Receiver';
[status,cmdout] = system(command,'-echo');
disp('Done -->> File is recovered');

%------------------- END ----------------------%
%----------------------------------------------%
