% This code simulates the three cases described in section 2.8.4 and
% produces the probability of interference versus distance graph in Figure 5
clear 
close all
clc

%% system model parameters
pf=.1; % target probability of false alarm 
pd=.9; % target probability of detection
R=500:2:2000; % distance between PU TX and each SU TX
iter=10000;
d0=1; % reference distance
n=3; % path loss exponent

R_pu=500; % primary user transmission radius
R_su=384.9326; % secondary user transmission radius

protection_factor=.95; 

p_n_dbm=-100;% noise power in dbm
p_n_db=-100-30;
p_n=10^((p_n_dbm-30)/10);% noise power in watt
snr_su_lim_db=3; % snr threhold for successful reception
snr_su_lim=10^(snr_su_lim_db/10);
snr_pu_lim_db=3;
snr_pu_lim=10^(snr_pu_lim_db/10);

p_pu=snr_pu_lim*p_n*(R_pu/d0)^n; % PU transmission power 
p_pu_db=10*log10(p_pu);

R_p=R_pu*(protection_factor); % protected radius

p_su=snr_su_lim*p_n*(R_su/d0)^n; % SU transmission power 
p_su_db=10*log10(p_su);

M = 4;                  % modulation index for psk
hpsk = comm.PSKModulator('ModulationOrder',M,...
    'BitInput',false,...
    'PhaseOffset',0);   % M-psk modulator


%% Case 1: Traditional CR
p_int_max=(p_pu*(R_p/d0).^-n/snr_pu_lim-p_n); % max interference power
% energy detector parameters
mu=(R_p/R_pu)^-n;
R_n=d0*((mu-1)*p_n/p_su)^-(1/n)+R_p; % no-talk radius
R_i=R_n-R_p;
snr_detect=0.1; % minimum detectable snr
snr_detect_db=10*log10(snr_detect);

N=ceil((qfuncinv(pf)-qfuncinv(pd)*sqrt(2*snr_detect+1)).^2./snr_detect.^2); % minimum samples required to perform energy detection

eps= (qfuncinv(pf)/sqrt(N)+1)*p_n; % detector threshold
eps_db=10*log10(eps);

R_su_pu_rx=abs(R-R_p); % distance between su and pu

p_pu_rx=p_pu*(R/d0).^(-n); 
snr_vect=p_pu_rx/p_n;
snr_vect_db=10*log10(snr_vect);

for i_snr=1:length(snr_vect_db)
    i_snr
    %% creating and modulating signal at PU transmitter
    infoSignal = randi(M,N,1)-1;  % random binary signal  (bits = log2(M))
    txSignal = sqrt(p_pu_rx(i_snr))*step(hpsk,infoSignal);   % M-psk signal
    temp_int=zeros(1,iter);
    gamma_n_temp=zeros(1,iter);
    for j = 1:iter     % no. of simulations
        j;
        rxSignal = awgn(txSignal,snr_vect_db(i_snr),'measured'); % AWGN channel
        
        p_r(j)=sum((abs(rxSignal)).^2)/(length(rxSignal)-1);%recieved energy
        
        if p_r(j) > eps % if energy is greater than threshold then signal is p_resent
            gamma_n_temp(j)=p_pu*(R_p/d0)^(-n)./(p_n);            
        else  temp_int(j)=p_su*(R_su_pu_rx(i_snr)/d0).^(-1*n);
            gamma_n_temp(j)=p_pu*(R_p/d0)^(-n)./(temp_int(j)+p_n);            
        end
    end
    prb_SINR(i_snr)=sum(gamma_n_temp>10^0.3)/iter;
end

%% Case 2: JT with traditional parameters

N_bs=2;
% energy detector parameters
OL_fac=.5; % amount of overlap between the base stations
d_bs=2*R_su-R_su*OL_fac; %% distance betweeen CBSs

R_pair=sqrt(R.^2-(d_bs/2)^2); % distance between PU TX and SU TX pair
R_su_pu_rx_pair=abs(R_pair-R_p);% distance between PU RX and SU TX pair
R_su_pu_rx_2bs=sqrt(R_su_pu_rx_pair.^2+(d_bs/2)^2); % distance between PU RX and each SU TX


M = 4;                  % modulation index for psk
hpsk = comm.PSKModulator('ModulationOrder',M,...
    'BitInput',false,...
    'PhaseOffset',0);   % M-psk modulator

for i_snr=1:length(snr_vect_db)
    i_snr
    %% creating and modulating signal at PU transmitter
    infoSignal = randi(M,N,1)-1;  % random binary signal  (bits = log2(M))
    txSignal = sqrt(p_pu_rx(i_snr))*step(hpsk,infoSignal);   % M-psk signal
    
    %% Additional codes for system model and interference power
    temp_int=zeros(1,iter);
    gamma_n_temp=zeros(1,iter);
    for j = 1:iter     % no. of simulations
        j;
        for i=1:N_bs
            rxSignal = awgn(txSignal,snr_vect_db(i_snr),'measured'); % AWGN channel            
            p_r=sum((abs(rxSignal)).^2)/(length(rxSignal)-1);%recieved energy            
            if p_r <eps % if energy is greater than threshold then signal is present                
                temp_int(j)=temp_int(j)+p_su*(R_su_pu_rx_2bs(i_snr)/d0).^(-1*n);             
            end 
        end
    end
    gamma_n_temp=p_pu*(R_p/d0)^(-n)./(temp_int+p_n);
    prb_SINR_2bs_traditional(i_snr)=sum(gamma_n_temp>10^0.3)/iter; 
end

%% Case 3: JT with CoMP parameters
R_i_2bs=d0*((mu-1)*p_n/2/p_su)^(-1/n);

x=sqrt(R_i_2bs^2-(d_bs/2)^2);
R_n_2bs=sqrt(R_i_2bs^2+R_p^2+2*R_p*x);
snr_detect_2bs=p_pu*(R_n_2bs/d0)^(-n)/p_n;
snr_detect_2bs_db=10*log10(snr_detect_2bs);
N_2bs=ceil((qfuncinv(pf)-qfuncinv(pd)*sqrt(2*snr_detect_2bs+1)).^2./snr_detect_2bs.^2);

eps= (qfuncinv(pf)/sqrt(N_2bs)+1)*p_n;
eps_db=10*log10(eps);

for i_snr=1:length(snr_vect_db)
    i_snr
    %% creating and modulating signal at PU transmitter
    infoSignal = randi(M,N_2bs,1)-1;  % random binary signal  (bits = log2(M))
    txSignal = sqrt(p_pu_rx(i_snr))*step(hpsk,infoSignal);   % M-psk signal
    
    
    temp_int=zeros(1,iter);
    gamma_n_temp=zeros(1,iter);

    for j = 1:iter     % no. of simulations
        j;
         for i=1:N_bs
            rxSignal = awgn(txSignal,snr_vect_db(i_snr),'measured'); % AWGN channel            
            p_r(j)=sum((abs(rxSignal)).^2)/(length(rxSignal)-1);%recieved energy
            if p_r(j) <eps % if energy is greater than threshold then signal is p_resent             
                temp_int(j)=temp_int(j)+p_su*(R_su_pu_rx_2bs(i_snr)/d0).^(-1*n); 
            end           
        end  
        gamma_n_temp(j)=p_pu*(R_p/d0)^(-n)./(temp_int(j)+p_n); 
    end
    prb_SINR_2bs_comp(i_snr)=sum(gamma_n_temp>10^0.3)/iter;
end

%% plotting probability vs distance
figure(1)
hold on
sk=5;

xlabel('Distance between SU TX and PU RX(m) ','FontSize',12)
ylabel('P( SINR > \gamma_{dec}  )','FontSize',12)

plot(R_su_pu_rx(1:sk:end),prb_SINR(1:sk:end), 'ks','linewidth', 2)
plot(R_su_pu_rx_2bs(1:sk:end),prb_SINR_2bs_traditional(1:sk:end),'r^','linewidth', 2)
plot(R_su_pu_rx_2bs(1:sk:end),prb_SINR_2bs_comp(1:sk:end), 'bo','linewidth', 2)

plot(R_su_pu_rx,prb_SINR, 'k','linewidth', 2)
plot(R_su_pu_rx_2bs,prb_SINR_2bs_traditional,'r','linewidth', 2)
plot(R_su_pu_rx_2bs,prb_SINR_2bs_comp, 'b','linewidth', 2)

plot(R_su_pu_rx,pd*ones(1,length(R_su_pu_rx)),'r--')

figure(1)
legend(['Traditional CR                ';
        'JT with Traditional parameters';
        'JT with CoMP parameters       '])
ylim([0.85 1.01]) 
xlim([R_su_pu_rx(1) R_su_pu_rx_2bs(end)])
grid on
box on    

