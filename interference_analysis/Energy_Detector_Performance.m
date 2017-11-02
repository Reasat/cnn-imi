% This code determine energy detector threshold and minimum number of 
% samples required to achieve a target probability of false detection and
% probability of false alarm (Section 2.8.3). Then it tests the CoMP-JT 
% energy detector performance at various SNRs and produce Figure 4
clear all
close all 
clc

snr_db= -11.9384 ; % minimun detectable snr
snr=10^(snr_db/10);
pf=.1; % target probability of false alarm 
pd=.9; % target probability of detection

p_n=1e-13; % noise power
p_n_db=10*log10(p_n);
p_n_dbm=p_n_db+30;

N=ceil((qfuncinv(pf)-qfuncinv(pd)*sqrt(2*snr+1)).^2./snr.^2); % minimum number of samples   
eps=(qfuncinv(pf)/sqrt(N)+1)*p_n; % detection threshold
eps_db=10*log10(eps);

snr_db_vect=-15:1:5;

p_s_db=snr_db_vect+p_n_db;
p_s_dbm=p_s_db+30;%signal power in dbm
p_s=10.^((p_s_dbm-30)/10); % power in watt/linear scale
snr_vect=p_s/p_n;

M = 4;                  % modulation index for psk
    hpsk = comm.PSKModulator('ModulationOrder',M,...
        'BitInput',false,...
        'PhaseOffset',0);   % M-psk modulator
iter=10000;    

for i_snr=1:length(snr_db_vect)
    i_snr
    %% creating and modulating signal at PU transmitter  
    infoSignal = randi(M,N,1)-1;  % random binary signal  (bits = log2(M))
    txSignal = sqrt(p_s(i_snr))*step(hpsk,infoSignal);   % M-psk signal
         
    %% determining pd experimentally
    
    d=0;
    for j = 1:iter     % no. of simulations
        j;
        rxSignal = awgn(txSignal,snr_db_vect(i_snr),'measured'); % AWGN channel
        p_r(j)=sum((abs(rxSignal)).^2)/(length(rxSignal)-1);%recieved energy
        if p_r(j) > eps % if energy is greater than threshold then signal is p_resent
            d = d+1;
        end
    end
    
    pd_ex(i_snr) = d/iter; % avg over 1000 simulation
    %% determining pd theoretically
    pd_th(i_snr)=qfunc((eps/p_n-snr_vect(i_snr)-1)*sqrt(N/(2*snr_vect(i_snr)+1)));
    %% determining pf experimentally
    d=0;
    for j = 1:iter     % no. of simulations
        j;
        n=(randn(1,N)+1i*randn(1,N))./(sqrt(2)); % Primary User Gaussian Signal
        rxSignal = sqrt(p_n)*n; % only noise, noise power in dbW
        p_r(j)=sum((abs(rxSignal)).^2)/(length(rxSignal)-1);%recieved energy
        if p_r(j) > eps % if energy is greater than threshold then signal is p_resent
            d = d+1;
        end
    end
    pf_ex(i_snr) = d/iter; % avg over 1000 simulation
    
    
    %% determining pf theoretically
    pf_th(i_snr) = qfunc((eps/p_n-1)*sqrt(N));
end

figure(1)
hold on

sk1=2;
sk2=4;
h1=plot(snr_db_vect(1:sk1:end),pd_ex(1:sk1:end),'bs','linewidth',2);
h2=plot(snr_db_vect(1:sk1:end),pf_ex(1:sk1:end),'rd','linewidth',2);
h3=plot(snr_db_vect(1:sk2:end),pd_th(1:sk2:end),'^','color',[132 45 1]/255,'linewidth',2);
h4=plot(snr_db_vect(1:sk2:end),pf_th(1:sk2:end),'ko','linewidth',2);
plot(snr_db_vect,pd_ex,'b','linewidth',2)
plot(snr_db_vect,pf_ex,'r','linewidth',2)
plot(snr_db_vect,pd_th,'color',[132 45 1]/255,'linewidth',2)
plot(snr_db_vect,pf_th,'k','linewidth',2);
xlabel('SNR (dB)','FontSize',12)
ylabel('P_d and P_f','FontSize',12)
ylim([0 1.1])
legend([h1,h2,h3,h4],{'P_d (simulated)','P_f (simulated)','P_d (analytical)','P_f (analytical)'})
box on

