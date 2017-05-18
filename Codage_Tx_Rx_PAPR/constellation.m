function constellation(freq_cor_IQData,nbsym)

LTS= [0.0000 + 0.0000i
  -0.0000 + 0.0000i
   0.0000 + 0.0000i
  -0.0000 - 0.0000i
   0.0000 + 0.0000i
  -0.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 - 0.0000i
  -8.0000 + 0.0000i
  -8.0000 - 0.0000i
   8.0000 + 0.0000i
   8.0000 + 0.0000i
  -8.0000 - 0.0000i
   8.0000 - 0.0000i
  -8.0000 - 0.0000i
   8.0000 + 0.0000i
   8.0000 - 0.0000i
   8.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 - 0.0000i
  -8.0000 + 0.0000i
  -8.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 + 0.0000i
  -8.0000 - 0.0000i
   8.0000 + 0.0000i
  -8.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 + 0.0000i
   8.0000 - 0.0000i
   8.0000 + 0.0000i
   0.0000 - 0.0000i
   8.0000 - 0.0000i
  -8.0000 + 0.0000i
  -8.0000 - 0.0000i
   8.0000 - 0.0000i
   8.0000 - 0.0000i
  -8.0000 + 0.0000i
   8.0000 - 0.0000i
  -8.0000 - 0.0000i
   8.0000 + 0.0000i
  -8.0000 - 0.0000i
  -8.0000 - 0.0000i
  -8.0000 + 0.0000i
  -8.0000 + 0.0000i
  -8.0000 + 0.0000i
   8.0000 - 0.0000i
   8.0000 + 0.0000i
  -8.0000 - 0.0000i
  -8.0000 - 0.0000i
   8.0000 - 0.0000i
  -8.0000 + 0.0000i
   8.0000 - 0.0000i
  -8.0000 - 0.0000i
   8.0000 + 0.0000i
   8.0000 - 0.0000i
   8.0000 - 0.0000i
   8.0000 - 0.0000i
   0.0000 + 0.0000i
   0.0000 - 0.0000i
  -0.0000 - 0.0000i
  -0.0000 - 0.0000i
  -0.0000 + 0.0000i];

% Received long training sequence

rlts=freq_cor_IQData(161+32:320);
rlts1=rlts(1:64);
rlts2=rlts(65:128);

RLTS1=fft(rlts1);
RLTS1=[RLTS1(33:64); RLTS1(1:32)];

RLTS2=fft(rlts2);
RLTS2=[RLTS2(33:64); RLTS2(1:32)];

%CFO_est = angle(Y(2,:)*Y(1,:)')/(2*pi);


%LS channel estimation
Hh_LS1=RLTS1./LTS;
Hh_LS2=RLTS2./LTS;

null_carriers=[1 2 3 4 5 6 33 60 61 62 63 64];
Hh_LS1(null_carriers)=0;%RLTS1(null_carriers);
Hh_LS2(null_carriers)=0;%RLTS2(null_carriers);

Y2=RLTS2';
Y1=RLTS1';

CFO_est = (angle(Y2*Y1')/(2*pi))*(1/(64*0.05e-6))

Hh_LS=(Hh_LS1+Hh_LS2)/2;



% MMSE
nFFT = 64; 
No = 0;
for i=1:nFFT
    temp(i)=(conj(RLTS1(i)-RLTS2(i))*(RLTS1(i)-RLTS2(i)));
    No = No +temp(i);
end
No=No/(nFFT*2);

X=diag(LTS);
Y=(RLTS1+RLTS2)/2;

autocor=xcorr(Hh_LS);
%plot(real(autocor))
%Rhh=diag(autocor(1:64));
Rhh=toeplitz(autocor(1:64));
F = dftmtx(nFFT)/sqrt(nFFT);

Rhy = Rhh * F' * X';        %X' is Hermitian Conjugate and Transpose of X
Ryy = X * F * Rhh * F' * X' + No * eye(nFFT);
Hh_MMSE = F * Rhy * inv(Ryy) * Y;
Hh_MMSE(null_carriers)=0;%RLTS1(null_carriers);

% channel estimation method selection
Hhat=Hh_LS;
%Hhat=Hh_MMSE;

y=freq_cor_IQData;
y=y(161:80*(nbsym+4));                 % Removing the short training sequence adding 4 symbols to compensate for the 4 symbols removed (2 short training seq and 2 long training sequence)                  
LTS=y(1:160);                          % Long training sequence with CP
LP=LTS(33:160);                        % Long training sequence without CP
LP1=LP(1:64);                          % 1st symbol
LP2=LP(65:128);                        % 2nd symbol

yycp=reshape(y,80,nbsym+2);             % with CP
yy=zeros(64,nbsym+2);                   % without CP
yy(:,1)=LP1;
yy(:,2)=LP2;
yy(:,3:end)=yycp(17:end,3:end);       % without CP


YY=fft(yy);
YY(:,:)=[YY(33:64,:) ;YY(1:32,:)];

Xp=[8+0*j;8+0*j; 8+0*j; -8+0*j];
kk=[12 26 40 54];
Yp1 = YY(:,3);
Yp2 = YY(:,4);
CFO_est = angle(sum(Yp2(kk).*Xp.*conj(Yp1(kk).*Xp)))/(2*pi)*64/80

YY=YY(:,3:nbsym+2) ;                   % no need for the LTS symbols
YY=reshape(YY,64*(nbsym),1);


k=0;
YY_eq=[];
YY_p_eq=[];
for m=1:nbsym
    YY_eq1=YY(1+k:64+k)./Hhat;
    YY_eq1(null_carriers)=0;%YY(null_carriers);
    YY_p_eq1=YY_eq1(kk);
    YY_p_eq=[YY_p_eq;YY_p_eq1];
    YY_eq=[YY_eq ;YY_eq1];
    k=k+64;
end
%CFO_est = angle(Yp2(kk).*Xp*(Yp1(kk).*Xp)')/(2*pi)*64/80;

%=================================== Correction de phase =====================================
%Charger le fichier des pilotes utilis�s lors de la transmission
itload sent_pilots_file.it
YY_p_th=[1;1;1;-1;1;1;1;-1;
sent_pilots(1:22*4)]*8;

phase_p=[];
phase_pc=[];
phase_pt=[];
YY_eq_pp = YY_eq;
for m=0:nbsym-1
    Ype = YY_eq(m*64+1:m*64+1+63);
    Yp = YY_p_th(m*4+1:m*4+4);
%save th_pilot YY_p_eq;
%Ype1 = YY_eq(1*64+1:1*64+1+63)';
%Ype2 = YY_eq(2*64+1:2*64+1+63)';
%CFO_est = angle(Yp2(kk).*Xp*(Yp1(kk).*Xp)')/(2*pi)*64/80;
%Moyennage sur les 4 pilotes du symbole OFDM
phase_p1=angle(sum(Yp.*conj(Ype(kk))));
%Pas de moyennage
phase_pt1=angle(Yp.*conj(Ype(kk)));

phase_p =[phase_p phase_p1];
phase_pt =[phase_pt; phase_pt1];

%On applique la correction de phase
YY_eq_pp(m*64+1:m*64+1+63)=YY_eq(m*64+1:m*64+1+63)*exp(j*phase_p1);

%On v�rifie que cela marche...
YY_eq_test=YY_eq_pp(m*64+1:m*64+1+63);
phase_p1c=angle(sum(Yp.*conj(YY_eq_test(kk))));
phase_pc =[phase_pc phase_p1c];
%phase_t2=angle(Yp2(kk)*conj(Ype2(kk))')
%scatterplot(Ype/64)
%pause;
end
%Affichage phase moyenn�e par symbole OFDM
phase_pd=(phase_p*180)/pi;
%Affichage phase moyenn�e trame corrig�e
phase_pdc=(phase_pc*180)/pi;
figure(4),plot(phase_pd);
xlabel('Symbole OFDM');
ylabel('Phase (degr�s)');
grid
hold
plot(phase_pdc,'r')
hold off
%Affichage phase instantan�e non moyenn�e
figure(5),plot((phase_pt*180)/pi)
xlabel('Pilotes');
ylabel('Phase (degr�s)');
grid
%pause;
YY=YY/sqrt(64);     % to display them
YY_eq=YY_eq/sqrt(64);
YY_eq_pp=YY_eq_pp/sqrt(64);

figure(2)
subplot(2,1,1);scatter(real(YY),imag(YY),'.','b');
axis([-2 2 -2 2]);
grid on;xlabel ('I');ylabel ('Q');title ('Received Data constellation before equalization')
subplot(2,1,2);scatter(real(YY_eq),imag(YY_eq),'.','r');
axis([-1.3 1.3 -1.3 1.3]);
grid on;xlabel ('I');ylabel ('Q');title ('Received Data constellation after equalization')
figure(3)
subplot(2,1,1);scatter(real(YY),imag(YY),'.','b');
axis([-2 2 -2 2]);
grid on;xlabel ('I');ylabel ('Q');title ('Received Data constellation before equalization')
subplot(2,1,2);scatter(real(YY_eq_pp),imag(YY_eq_pp),'.','r');
axis([-1.3 1.3 -1.3 1.3]);
grid on;xlabel ('I');ylabel ('Q');title ('Received Data constellation after equalization and phase correction')


% Phase tracking
% Mehod is proposed by this paper "FPGA implementation of an OFDM-based
% WLAN receiver" by Mar?a Jos� Canet ?, Javier Valls, Vicen� Almenar, Jos� Mar?n-Roig
% 
% !!! DOESNOT WORK WELL
% Polarity_Sequence = [1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1,1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1,-1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
% Pilots = [1+j*0 1+j*0 1+j*0 -1+j*0];
% pilot_index=[11 25 39 53]+1;  % 1  : matlab  index starts from 1
% k=0;                          % 0 for the signal field
% for m=1:64:64*(nbsym)
%    
%     for ind=1:4
%           P(ind)=(Pilots(ind)*Polarity_Sequence(mod(k,126)+1));
%           R(ind)=YY_eq(pilot_index(ind)+m-1);
%           Hp(ind)=Hh_LS(pilot_index(ind));
%     end
%     P
%     R
%     Hp
%     conj(Hp.*P)
%     R.*conj(Hp.*P)
%     sum(R.*conj(Hp.*P))
%     angle(sum(R.*conj(Hp.*P)))
%   if sign(real(P))~=sign(real(R))
%     k
%     P
%     R
%  
%   end
%   n=k+1;
%     mult(n)=sum(R.*conj(Hp.*P));   
%     phase(n)=angle(mult(n));
%    
%     if n>=2
%         phase(n-1)-phase(n);
%     end
%     
%     k=k+1;
% 
% end  
% phase;
% YY_track=[];
% tot_phase=0;
% for k=1:(nbsym)
%     smybol=1+(k-1)*64:k*64;
%     phase(k)/(2.65e4)*smybol
%    YY_track(smybol)=YY_eq(smybol).*exp((-j*phase(k)/(2.65e4))*smybol)';
%     tot_phase=tot_phase+phase(k);
%    YY_track(smybol)=YY_eq(smybol)*exp(-j*tot_phase/(4.18*1e2));
%     tot_phase
% end
% 
% YY_track=YY_track*exp(-j*phase(1));
% 
% subplot(3,1,3);scatter(real(YY_track),imag(YY_track),'.','b');
% grid on;xlabel ('I');ylabel ('Q');title ('Received Data constellation after phase tracking')
% axis([-1.2 1.2 -1.2 1.2]);
% 
% 
% itload Sent_Frame.it;           % loading the IQ file
% IQData=PPDU_Frame_Final(1:end);
% 
% 
% x=IQData;
% x=x(161:end);                         
% nbsym=length(x)/80;
% LP=x(1:160);                          % LP with CP
% LP=x(33:160);                         % LP without CP
% LP1=LP(1:64);                         % 1st Long training sequence
% LP2=LP(65:128);                       % 2nd Long training sequence
% 
% xxcp=reshape(x,80,nbsym);             % with CP
% xx=zeros(64,nbsym);                   % without CP
% xx(:,3:end)=xxcp(17:end,3:end);
% xx(:,1)=LP1;
% xx(:,2)=LP2;
% XX=fft(xx);
% XX(:,:)=[XX(33:64,:) ;XX(1:32,:)];
% XX=XX(:,3:nbsym) ;          % no need for the LTS symbols
% XX=reshape(XX,64*(nbsym-2),1);
% XX=XX/sqrt(64);

% figure(3)
% for i=1:1000 % randomly
%     scatter(real(XX(i)),imag(XX(i)),'.','r');hold on 
%     scatter(real(YY_eq(i)),imag(YY_eq(i)),'.','b')
%     axis([-1.2 1.2 -1.2 1.2]);grid on;xlabel ('I');ylabel ('Q');
%     pause(.05)
%     hold off;
% end
