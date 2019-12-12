function [tfr_noise,tfr_noise_int,Band_noise] = compute_tfr(sig,window,SNR)
 
 %INPUT
 %sig        : type of studied signal
 %window     : window used 
 %SNR        : SNR corresponding to the added noise
 
 %OUTPUT
 %tfr_noise     : noisy STFT
 %tfr_noise_int : thresholded STFT
 %Band_noise    : Band of coefficient kept for noise interpolation
 
 N = 4096;
 if sig == 1
  t  = (0:N-1)/N;
  a  = 2;
  s1 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
  s2 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
  s  = s1+s2;
  s = s(:);
  Sacc = [s1(:) s2(:)];
  nr = 2;
  Lg = 161;
  sigma_opt = 0.15;
  clwin = 10;
 elseif sig == 2 
  s  = fmlin(N,.05,0.25);
  Sacc = s;
  nr = 1;
  Lg = 161;
  sigma_opt = 0.15;
  clwin = 10;
 end
 
 if (sig <= 2)
  Nfft = 512;
 end
 
 s    = s(:);

 %we build the filter h
 if strcmp(window,'hamming')
  hlength=floor(Lg);%optimal window determined by Rényi entropy
  hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
  h = tftb_window(hlength,window);
  [hrow,hcol]=size(h); 
  Lh=(hrow-1)/2;
 else
  %the window is the Gaussian window    
  prec = 10^(-3);
  L =  sigma_opt*Nfft;
  Lh = floor(L*sqrt(-log(prec)/pi))+1;
  h = amgauss(2*Lh+1,Lh+1,L); 
 end
 %tfr sans bruit
 [tfr0,norm2h] = tfrstft_three_case_down(s,Nfft,2,h,Lh,1,0); 

 n = randn(N,1)+1i*randn(N,1);
 [sn] = sigmerge(s,n,SNR);
 sn   = sn(:);
 
 %tfr debruitée globale
 [tfr_noise,norm2h] = tfrstft_three_case_down(sn,Nfft,2,h,Lh,1,0);
 Y2 = real(tfr_noise);
 gamma_estime = median(abs(Y2(:)))/0.6745;
 [Cs] = exridge_mult(tfr_noise,nr,0,0,clwin);
 Cs = Cs';
 
 B = size(tfr_noise); 
 interval  = zeros(B(2),nr);
 Abstfr = abs(tfr_noise);
 
 for j=1:nr,
   %construction of the TF mask
   tfr_noise_int = zeros(B); 
   Band_noise    = ones(B);
   for r = 1:B(2),  
   
   val = 3*gamma_estime; %trheshold for the transfrom depending on the noise level
  
   if (Abstfr(Cs(r,j),r) > val)
    k1 = 0;
    k2 = 0;
     
    eta1 = - 1;
    while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > val)
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
    while (eta2 < 0) && (Abstfr(Cs(r,j)+min(k2,B(1)-Cs(r,j)),r) > val)
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
    tfr_noise_int(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r) = tfr_noise(max(1,Cs(r,j)-eta1):min(B(1),Cs(r,j)+eta2),r);
    Band_noise(max(1,Cs(r,j)-eta1-q):min(B(1),Cs(r,j)+eta2+q),r) = 0;
   end
  end
 end
 
   
 

 
 
  
 