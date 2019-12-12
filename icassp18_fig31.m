function icassp18_fig31()
% plots figure 3a of paper "NOVEL THRESHOLDING TECHNIQUE FOR THE DENOISING OF MULTICOMPONENT
%SIGNALS" by Duong-Hung and Meignen.

clear all; close all; clc; 
set(0,'DefaultAxesFontSize',20);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/ICASSP2018/figures';
N = 4096;

t = (0:N-1)';
Nfft = 512;
freq = (0:Nfft-1)/Nfft;

[tfr_freenoise,tfr_noise,tfr_noise_hard,tfr_noise_soft,gamma_estime,s,h,Lh,sn] = compute_tfr_bis(1,'Gauss',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = size(tfr_freenoise);


FigHandle(1) = figure(); colormap(1-gray);
imagesc(t,freq,abs(tfr_freenoise(1:B(1)/2,:))); %xlim([0.1 0.9]);
set(gca,'ydir','normal');% title('noisy STFT') 
 ylabel('\omega');xlabel('L'); %xlim([0.1*N 0.9*N])
axis square

pause

for i=1:1
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/icassp_fig21', chemin0),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end
end

