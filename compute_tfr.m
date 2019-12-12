function [tfr_freenoise,tfr_noise,tfr_noise_hard,tfr_noise_soft,s,h,Lh,sn,gamma_estime] = compute_tfr1(sig,window,SNR)

%INPUT
%sig        : type of studied signal
%window     : window used 
%SNR        : SNR corresponding to the added noise

%OUTPUT
%tfr_noise_hard     : noisy STFT, hard thresholded
%tfr_noise_soft     : noisy STFT, soft thresholded

N = 4096;
if sig == 1
    t  = (0:N-1)/N;
    a  = 2;
    s1 = a.*exp(2*pi*1i*(600*t+30*cos(3*pi*t)));
    s2 = a.*exp(2*pi*1i*(1200*t+60*cos(3*pi*t)));
    s  = s1+s2;
    s = s(:);
    Sacc = [s1(:) s2(:)];
    nr = 2;
    Lg = 161;
    sigma_opt = 0.15;
    clwin = 10;
elseif sig == 2 
    s  = fmlin(N,.1,0.35);
    %s = fmconst(N,0.5);
    Sacc = s;
    nr = 1;
    Lg = 161;
    sigma_opt = 0.15;
    clwin = 10;
elseif sig == 3 
    t  = (0:N-1)/N;
    a  = 2;
    s = a.*exp(2*pi*1i*(1500*t+30*cos(3*pi*t)));
    s = s(:);
    Sacc = s;
    nr = 1;
    Lg = 161;
    sigma_opt = 0.15;
    clwin = 10;
elseif sig == 4 
    load -ascii batsig.txt
    s = [batsig];
    s=hilbert(s);
    N = length(s);
    nr = 2;
    sigma_opt = 0.15;
    clwin = 10;
end

if (sig <= 3)
    Nfft = 512;
else
    Nfft = N/2;   
end
s    = s(:);
%we build the filter h
if strcmp(window,'hamming')
    hlength=floor(Lg);%optimal window determined by Renyi entropy
    hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
    h = tftb_window(hlength,window);
    [hrow,hcol]=size(h); 
    Lh=(hrow-1)/2;
else
    %the window is the Gaussian window    
    prec = 10^(-6);
    L =  sigma_opt*Nfft;
    Lh = floor(L*sqrt(-log(prec)/pi))+1;
    h = amgauss(2*Lh+1,Lh+1,L); 
end

n = randn(N,1)+1i*randn(N,1);
[sn] = sigmerge(s,n,SNR);
sn   = sn(:);

[tfr_noise,norm2h] = tfrstft_three_case_down(sn,Nfft,2,h,Lh,1,0);
[tfr_freenoise,~] = tfrstft_three_case_down(s,Nfft,2,h,Lh,1,0);

Y2 = real(tfr_noise);
gamma_estime = median(abs(Y2(:)))/0.6745;

[Cs] = exridge_mult(tfr_noise,nr,0,0,clwin);

Cs = Cs';

B = size(tfr_noise); 
interval  = zeros(B(2),nr);
Abstfr      = abs(tfr_noise);

tfr_noise_shift_real = zeros(B);
tfr_noise_shift_imag = zeros(B);
tfr_noise_shift  = zeros(B);

tfr_noise_sym= zeros(B);

for j=1:nr
    %construction of the TF mask
    tfr_noise_hard = zeros(B); 
    tfr_noise_soft = zeros(B);
    for r = 1:B(2)
        val = 3*gamma_estime; %threshold for the transform depending on the noise level
        if (Abstfr(Cs(r,j),r) > sqrt(2)*val)
            k1 = 0;
            k2 = 0;
            eta1 = - 1;
            while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) >sqrt(2)*val)
                if (k1 ~= Cs(r,j)-1)
                    k1 = k1+1;
                else
                    eta1 = k1;
                end
            end
            if (eta1 < 0)
                eta1 = k1-1;
            end
            
            eta2 = -1;
            while (eta2 < 0) && (Abstfr(Cs(r,j)+min(k2,B(1)-Cs(r,j)),r) > sqrt(2)*val)
                if (k2 ~= B(1)-Cs(r,j))
                    k2 = k2+1;   
                else
                    eta2 = k2;
                end
            end
            if (eta2 < 0)
                eta2 = k2;
            end
            interval(r,j) =eta2+eta1+1;
            eta = max(eta1,eta2); %we take the larger value
            
            indmax = max(1,Cs(r,j)-eta);
            indmin  = min(B(1),Cs(r,j)+eta);
            
            tfr_noise_hard(indmax:indmin,r) = tfr_noise(indmax:indmin,r); 
            
            [~,I] = max(abs(tfr_noise_hard(indmax:indmin,r)));
            [~,Ir] = max(abs(real(tfr_noise_hard(indmax:indmin,r))));
            [~,Ii] = max(abs(imag(tfr_noise_hard(indmax:indmin,r))));
            
%             figure(); 
%             hold on;
%             plot(-eta:eta,abs(tfr_noise_hard(indmax:indmin,r)),'o-','color','b')
%             plot(-eta:eta,real(tfr_noise_hard(indmax:indmin,r)),'o--','color','b')
%             plot(-eta:eta,imag(tfr_noise_hard(indmax:indmin,r)),'.-','color','b')
%             plot(-eta:eta,abs(tfr_freenoise(indmax:indmin,r)),'-+','color','r')
%             plot(-eta:eta,real(tfr_freenoise(indmax:indmin,r)),'--+','color','r')
%             plot(-eta:eta,imag(tfr_freenoise(indmax:indmin,r)),'.-','color','r')
%            
            tfr_noise_shift_real(:,r) = circshift(real(tfr_noise_hard(:,r)),sign(Ir-I)*(I-Ir));
            tfr_noise_shift_imag(:,r) = circshift(imag(tfr_noise_hard(:,r)),sign(Ii-I)*(I-Ii));            
            tfr_noise_shift(indmax:indmin,r) = tfr_noise_shift_real(indmax:indmin,r)+1i*tfr_noise_shift_imag(indmax:indmin,r);

                  
%             figure(); 
%             hold on;
%             plot(-eta:eta,abs(tfr_noise_shift(indmax:indmin,r)),'o--','color','b')
%             plot(-eta:eta,real(tfr_noise_shift(indmax:indmin,r)),'x-','color','b')
%             plot(-eta:eta,imag(tfr_noise_shift(indmax:indmin,r)),'-','color','b')
%             plot(-eta:eta,abs(tfr_freenoise(indmax:indmin,r)),'-+','color','r')
%             plot(-eta:eta,real(tfr_freenoise(indmax:indmin,r)),'--+','color','r')
%             plot(-eta:eta,imag(tfr_freenoise(indmax:indmin,r)),'.-','color','r')
            
            tfr_noise_sym(indmax:indmin,r) = 1/2*(tfr_noise_shift(indmax:indmin,r)+tfr_noise_shift(indmin:-1:indmax,r)); 
          
%             figure(); 
%             hold on;
%             plot(-eta:eta,abs(tfr_noise_sym(indmax:indmin,r)),'o--','color','b')
%             plot(-eta:eta,real(tfr_noise_sym(indmax:indmin,r)),'x-','color','b')
%             plot(-eta:eta,imag(tfr_noise_sym(indmax:indmin,r)),'-','color','b')
%             plot(-eta:eta,abs(tfr_freenoise(indmax:indmin,r)),'-+','color','r')
%             plot(-eta:eta,real(tfr_freenoise(indmax:indmin,r)),'--+','color','r')
%             plot(-eta:eta,imag(tfr_freenoise(indmax:indmin,r)),'.-','color','r')
            
            X  = [1:max(1,Cs(r,j)-eta-2) indmax:indmin min(B(1),Cs(r,j)+eta+2):B(1)];
            X  = unique(X);
            Y  = tfr_noise_sym(X,r);
            YY = pchip(X,real(Y),1:B(1));
            ZZ = pchip(X,imag(Y),1:B(1));
            tfr_noise_soft(:,r) = YY' + 1i*ZZ';                            
        end
    end
end
end
 
  
 