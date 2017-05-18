% Test 

itload sent_pilots_file.it

itload Sent_Frame.it

YY_p_th=[1;1;1;-1;1;1;1;-1;sent_pilots(1:22*4)']*8;

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
