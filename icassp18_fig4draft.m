% plots figure 4 of paper "NOVEL THRESHOLDING TECHNIQUE FOR THE DENOISING OF MULTICOMPONENT
%SIGNALS" by Duong-Hung and Meignen.

clear all; close all; clc; 
set(0,'DefaultAxesFontSize',20);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/ICASSP2018/figures';
myfilenameSNR = sprintf('SNRGain_bp_1.mat');
load SNRGain_BT.mat
SNRgain0 = importdata(myfilenameSNR);
SNRgainh = SNRgain0.SNRoutputGainh; SNRgainh = SNRgainh(:);
SNRgains = SNRgain0.SNRoutputGains; SNRgains = SNRgains(:);

 FigHandle(1) = figure(); %set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
 boxplot([SNRgain(:),SNRgainh,SNRgains],'Labels',{'BT','HT','SSR-HT'})
 ylabel('SNR Gain [dB]')

%   FigHandle(1) = figure(); %set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
%  boxplot([SNRgain(:,1),SNRgain(:,3),SNRgain(:,3)],'Labels',{'BT','HT','IHT'})
%  ylabel('SNR Gain [dB]')
 for i=1
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/icassp_fig4_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
 end
 
% FigHandle(1) = figure();
% plot(SNR,SNRout2,'k',SNR,SNRout1,'r.--',SNR,SNRoutBT,'b--+');
% legend('IHT','HT','BT','Location','best');
% xlabel('SNR in (dB)');ylabel('SNR out (dB)'); xlim([-5 30])
% explot();