% Test for PGC signal of HT, IHT, and BT
clear all; close all; %clc; 
set(0,'DefaultAxesFontSize',16);
addpath(genpath('/users/these/phamdu/Dropbox/Working/LettreIEEE'));
addpath(genpath('/users/these/phamdu/Dropbox/Working/'))
addpath(genpath('/media/phamdu/USB/PCG_data/SISECFull'))
%addpath(genpath('/users/these/phamdu/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO2018_Imroved_Hard_Thresholding/programmes'))

%chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO2018_Imroved_Hard_Thresholding/figures';
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/ICASSP2018/figures';

%% PCG signal 
SNRoutputGainBT = zeros(16,1);
SNRoutputGainh = zeros(16,1);
SNRoutputGains = zeros(16,1);

cof = 1;
for m=1:length(cof)
    COF=cof(m)
    for k= 1:16
        X = sprintf('PCG is %dth:', k);
        disp(X)

        myfilenameclean = sprintf('S%d_Clean.mat', k);
        myfilenamenoise = sprintf('S%d.mat', k);
        myfilenameECG = sprintf('S%d_ECG.mat', k);

        sclean = importdata(myfilenameclean);
        snoise = importdata(myfilenamenoise);
        sECG = importdata(myfilenameECG);

        s_clean = sclean.PCG;
        s_noise = snoise.x;
        fs = snoise.fs;
        s_ECG = sECG;

        N1 =length(s_noise);

        s_clean = s_clean(1:N1);
        s_noise =s_noise(1:N1);
        s_ECG = s_ECG(1:N1);

        window = 'gauss';
        downsamp = 1; 
        nr = 1;
        %cof = 40
        index =  (200:N1-200)';

        SNRinput(k)= snr(s_clean(index),s_clean(index)-s_noise(index))

        %%
        [tfr_freenoise,tfr_noise,tfr_noise_hard_nr,tfr_noise_soft_nr,modes_hard,modes_soft,gamma_estime,s,h,Lh,Cs] = compute_tfr_bis2(s_noise,s_clean,window,nr,cof(m));
        % Blockthresholding
        % noisy signal
        sigma_noise = exp(-SNRinput(k)/20); % sigma
        %sn = real(s)+sig*randn(size(s));        
        [srecBT, ~, ~] = BlockThresholding(s_noise, 127, fs, sigma_noise); srecBT = srecBT(:); % BlockThresholding(fn, time_win, f_sampling, sigma_noise)
        SNRoutputBT(k) = snr(s_clean(index),s_clean(index)-srecBT(index))
        SNRoutputh(k)= snr(s_clean(index),s_clean(index)-modes_hard(index));
        SNRoutputs(k) = snr(s_clean(index),s_clean(index)-modes_soft(index));
        
        SNRoutputGainBT(k) = SNRoutputBT(k) - SNRinput(k)
        SNRoutputGainh(k) = SNRoutputh(k) - SNRinput(k)
        SNRoutputGains(k) = SNRoutputs(k) - SNRinput(k)
        clearvars -except SNRoutputGainBT SNRoutputGainh SNRoutputGains chemin0 cof(m) m cof
    end

    FigHandle(1) = figure(); %set(FigHandle(1),'units','normalized','outerposition',[0 0 1 1]);
    boxplot([SNRoutputGainBT,SNRoutputGainh,SNRoutputGains],'Labels',{'BT','HT','IHT'})
    ylabel('SNR Gain [dB]')

    %%
    for i=1
        export_fig(FigHandle(i), ... % figure handle
            sprintf('%s/icassp_fig41_%d', chemin0,i),... % name of output file without extension
            '-painters', ...      % renderer
            '-transparent', ...   % renderer
            '-pdf', ...         % file format
            '-r5000' );             % resolution in dpi
    end
end

