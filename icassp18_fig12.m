%function icassp18_fig1()
%plots figures 1-2 of paper "NOVEL THRESHOLDING TECHNIQUE FOR THE DENOISING OF MULTICOMPONENT
%SIGNALS" by Duong-Hung PHAM and Meignen.

clear all; close all; clc; 
set(0,'DefaultAxesFontSize',20);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/ICASSP2018/figures';
N = 4096;
dt = 1/N;
t = (0:N-1)';
Nfft = 512;
freq = (0:Nfft-1)/Nfft;
index = 100:140;

sig=2; SNR=0; window = 'Gauss';

[tfr_freenoise,tfr_noise,tfr_noise_hard,tfr_noise_soft,gamma_estime,s,h,Lh,sn] = compute_tfr_bis(2,'Gauss',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = size(tfr_noise);
 
FigHandle(1) = figure(); colormap(1-gray);
imagesc(t,freq,abs(tfr_noise(1:B(1)/2,:)));
set(gca,'ydir','normal');% title('noisy STFT')
line([N/2-1 0.5],[0 N/2-1],'color','red','Linewidth',2); ylabel('\omega');xlabel('L'); %xlim([0.1*N 0.9*N])
axis square

%%
j=N/2;
FigHandle(2) = figure(); %set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(2*index/Nfft,abs(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,abs(tfr_noise(index,j)),'x-','color','b')
plot(2*index/Nfft,abs(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,3*sqrt(2)*gamma_estime*ones(1,41),'c--','LineWidth',1.5)
legend('noise-free','noisy','HT','HT level')
xlabel('\omega');ylabel('|V_f^g|');  
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([0 80]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off

FigHandle(3) = figure(); %set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(2*index/Nfft,real(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,real(tfr_noise(index,j)),'x-','color','b')
plot(2*index/Nfft,real(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,3*gamma_estime*ones(1,41),'c--','LineWidth',1.5)
legend('noise-free','noisy','HT','real HT level')
xlabel('\omega'); ylabel('Re(V_f^g)');  
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([0 80]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off


FigHandle(4) = figure(); %set(FigHandle(4),'units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(2*index/Nfft,imag(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,imag(tfr_noise(index,j)),'x-','color','b')
plot(2*index/Nfft,imag(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,3*gamma_estime*ones(1,41),'c--','LineWidth',1.5)
legend('noise-free','noisy','HT','imag HT level')
xlabel('\omega');ylabel('Im(V_f^g)');  
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([-15 40]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off


%%
FigHandle(5) = figure(); 
hold on;
plot(2*index/Nfft,abs(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,abs(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,abs(tfr_noise_soft(index,j)),'--','color','b','LineWidth',1.5)
legend('noise-free','HT','SSR-HT')
xlabel('\omega');ylabel('|V_f^g|');  
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([0 80]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off

FigHandle(6) = figure();
hold on;
plot(2*index/Nfft,real(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,real(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,real(tfr_noise_soft(index,j)),'--','color','b','LineWidth',1.5)
legend('noise-free','HT','SSR-HT')
xlabel('\omega');ylabel('Re(V_f^g)'); 
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([-15 80]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off


FigHandle(7) = figure(); 
hold on;
plot(2*index/Nfft,imag(tfr_freenoise(index,j)),'o--','color','k')
plot(2*index/Nfft,imag(tfr_noise_hard(index,j)),'-','color','r')
plot(2*index/Nfft,imag(tfr_noise_soft(index,j)),'--','color','b','LineWidth',1.5)
legend('noise-free','HT','SSR-HT')
xlabel('\omega');ylabel('Im(V_f^g)'); 
xlim([2*round(100/Nfft,2) 2*round(140/Nfft,2)]); 
ylim([-15 40]); 
set(gca,'XTick',2*round(100/Nfft,2):round(20/Nfft,2):2*round(140/Nfft,2))
hold off

%%
pause 
for i=1:4
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/icassp_fig1_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end

for i=5:7
    export_fig(FigHandle(i), ... % figure handle
        sprintf('%s/icassp_fig1_%d', chemin0,i),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r5000' );             % resolution in dpi
end


