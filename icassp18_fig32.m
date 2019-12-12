% plots figure 3b of paper "NOVEL THRESHOLDING TECHNIQUE FOR THE DENOISING OF MULTICOMPONENT
%SIGNALS" by Duong-Hung and Meignen.

clear all; close all; clc; 
set(0,'DefaultAxesFontSize',22);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/ICASSP2018/figures';

SNR = -5:5:30;
P = length(SNR);
rep = 10;
SNRout1 = zeros(1,P); % Old method 
SNRout2 = zeros(1,P); % New method
SNRoutBT = zeros(1,P); % Block thresholding
%SNRoutSP = zeros(1,P); % SP
%SNRoutOSG = zeros(1,P); % OSG

index = 100:4096-100;
index1 = 100:150;
lam = 0.008;
MAX_ITER = 50;

for k=1:length(SNR)
    for l=1:rep
        SNR(k)
        % Old and new method
        [tfr_freenoise,tfr_noise,tfr_noise_hard,tfr_noise_soft,gamma_estime,s,h,Lh,sn]= compute_tfr_bis(2,'Gauss',SNR(k));
        B = size(tfr_noise);
        [srec1] = itfrstft_three_case_down(tfr_noise_hard,2,B(2),h,0);
        [srec2] = itfrstft_three_case_down(tfr_noise_soft,2,B(2),h,0);
        
        % Blockthresholding
        % noisy signal
        sig = exp(-SNR(k)/20); % sigma
        %sn = real(s)+sig*randn(size(s));        
        [srecBT, ~, ~] = BlockThresholding(real(sn), 20, 4096, sig);
        
        %Denoise Simple Prior
        %srecSP =DenoiseSimplePrior(real(sn),lam,h,Lh,MAX_ITER);
        
%         OSG Technique
%         Nit = 25;                       % Nit : number of iterations
%         K = 5;                          % K : group size parameter
%         lambda = 0.68 * sig;          % lambda: regularization parameter for OGS   
%         srecOSG = ogshrink(real(sn), K, lambda, Nit);       % Run the OGS algorithm
% 
%         figure(); hold on; 
%         plot(index1,real(srec1(index1)),'--',index1,real(srec2(index1)),'+-',index1,srecBT(index1)',index1,srecSP(index1),'*-')
%         legend('Old method','New Method','BT','SP');
        
        SNRout1(k) = SNRout1(k)+ snr(s(index),s(index)-srec1(index));
        SNRout2(k) = SNRout2(k)+ snr(s(index),s(index)-srec2(index));
        SNRoutBT(k) = SNRoutBT(k)-10*log(norm(real(s(index))-srecBT(index)')/norm(real(s(index)))); 
%       SNRoutSP(k) = SNRoutSP(k)- 10*log(norm(real(s(index))-srecSP(index))/norm(real(s(index)))); 
%       SNRoutOSG(k) = SNRoutOSG(k)- 10*log(norm(real(s(index))-srecOSG(index))/norm(real(s(index))));     
    end
end
SNRout1=SNRout1/rep;
SNRout2=SNRout2/rep;
SNRoutBT=SNRoutBT/rep;
%SNRoutSP=SNRoutSP/rep;
%SNRoutOSG=SNRoutOSG/rep;

% FigHandle(1) = figure();
% plot(SNR,SNRout2,'k',SNR,SNRout1,'r.--',SNR,SNRoutBT,'b--+');
% legend('IHT','HT','BT','Location','best');
% xlabel('SNR in (dB)');ylabel('SNR out (dB)'); xlim([-5 30])
% explot();

FigHandle(1) = figure();
plot(SNR,SNRout(2,:),'k',SNR,SNRout(1,:),'r.--',SNR,SNRout(3,:),'b--+');
legend('SSR-HT','HT','BT','Location','best');
xlabel('SNR in [dB]');ylabel('SNR out [dB]'); xlim([-5 30])
explot();
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1
 %%%%%%%%%%%%%%%%%%%%%% print Figures
 export_fig(FigHandle(i), ... % figure handle
     sprintf('%s/icassp_fig22', chemin0),... 
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
end