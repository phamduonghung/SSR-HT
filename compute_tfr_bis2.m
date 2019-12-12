function [tfr,tfr_noise,tfr_noise_hard_nr,tfr_noise_soft_nr,modes_hard,modes_soft,gamma_estime,s,h,Lh,Cs] = compute_tfr_bis2(s,s_clean,window,nr,cof)

%%% Testing for real PCG signal
%INPUT
%sig        : type of studied signal
%window     : window used 
%SNR        : SNR corresponding to the added noise

%OUTPUT
%tfr_noise_hard     : noisy STFT, hard thresholded
%tfr_noise_soft     : noisy STFT, soft thresholded

%s=s_noise;
N = length(s);
Nfft = 512;
Lg = 64;
sigma_opt = 0.071;
clwin = 10;
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
    prec = 10^(-3);
    L =  sigma_opt*Nfft;
    Lh = floor(L*sqrt(-log(prec)/pi))+1;
    h = amgauss(2*Lh+1,Lh+1,L); 
end


[tfr_noise,~] = tfrstft_three_case_down(s,Nfft,2,h,Lh,1,0);
[tfr,~] = tfrstft_three_case_down(s_clean,Nfft,2,h,Lh,1,0);

Y2 = real(tfr_noise);
gamma_estime = median(abs(Y2(:)))/0.6745;

%tfr_noise1=tfr_noise;
tfr_noise(50:end,:)=0;

B = size(tfr_noise);

[Cs] = exridge_mult(tfr_noise(1:50,:),nr,0,0,clwin);

Cs = Cs';

B = size(tfr_noise); 
Abstfr      = abs(tfr_noise);

tfr_noise_soft_nr = zeros(B(1),B(2),nr);
tfr_noise_hard_nr = zeros(B(1),B(2),nr);
modes_hard = zeros(B(2),nr);
modes_soft = zeros(B(2),nr);

for j=1:nr
    %construction of the TF mask
    tfr_noise_hard = zeros(B); 
    tfr_noise_soft = zeros(B);

    for r = 1:B(2)  
        val = cof*gamma_estime; %threshold for the transform depending on the noise level
        if (Abstfr(Cs(r,j),r) > sqrt(2)*val)
            k1 = 0;
            k2 = 0;
            eta1 = - 1;
            while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > sqrt(2)*val)
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
            eta = max(eta1,eta2); %we take the larger value
            indmin = max(1,Cs(r,j)-eta);
            indmax = min(B(1),Cs(r,j)+eta);
            
            tfr_noise_hard(indmin:indmax,r) = tfr_noise(indmin:indmax,r);  % hard thresholding 

            %Shift step 
            
            if (abs(tfr_noise(Cs(r,j),r))-abs(tfr_noise(max(1,Cs(r,j)-1),r)) >= 2) &&...
                    (abs(tfr_noise(Cs(r,j),r))-abs(tfr_noise(Cs(r,j)+1,r)) >= 2)...
                    && (Cs(r,j)-eta >= 1) && (Cs(r,j)+eta <= B(1))
                [Y,I] = max(abs(imag(tfr_noise(indmin:indmax,r))));
                shift = I - eta -1;
                tfr_noise_shift_imag = circshift(imag(tfr_noise(:,r)),shift);
                [Y,I] = max(abs(real(tfr_noise(indmin:indmax,r))));
                shift = I - eta -1;
                tfr_noise_shift_real = circshift(real(tfr_noise(:,r)),shift);
                tfr_noise_shift = tfr_noise_shift_real + 1i*tfr_noise_shift_imag;
                
                %Symmetry step
                tfr_noise_soft(indmin:indmax,r) = 1/2*(tfr_noise_shift(indmin:indmax)+tfr_noise_shift(indmax:-1:indmin));
            else
                tfr_noise_shift = tfr_noise(:,r);
                
                %Symmetry step
                if abs(tfr_noise_shift(Cs(r,j)+1)) >= abs(tfr_noise_shift(max(1,Cs(r,j)-1)))
                    for k = 0:eta
                        if ((Cs(r,j)-k) >= 1)&&((Cs(r,j)+k+1) <= B(1))
                            tfr_noise_soft(Cs(r,j)-k,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k)+tfr_noise_shift(Cs(r,j)+k+1));
                            tfr_noise_soft(Cs(r,j)+k+1,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k)+tfr_noise_shift(Cs(r,j)+k+1));
                        end
                    end
                else
                    for k = 0:eta
                        if ((Cs(r,j)-k-1) >= 1)&&((Cs(r,j)+k) <= B(1))   
                            tfr_noise_soft(Cs(r,j)-k-1,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k-1)+tfr_noise_shift(Cs(r,j)+k));
                            tfr_noise_soft(Cs(r,j)+k,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k-1)+tfr_noise_shift(Cs(r,j)+k));
                        end
                    end
                end
            end            
         
            %smoothing step
            X = [1:max(1,Cs(r,j)-eta-4) max(1,Cs(r,j)-eta):min(B(1),Cs(r,j)+eta) min(B(1),Cs(r,j)+eta+4):B(1)];
            X  = unique(X);
            Y  = tfr_noise_soft(X,r);
            YY = pchip(X,real(Y),1:B(1));
            ZZ = pchip(X,imag(Y),1:B(1));
            tfr_noise_soft(:,r) = YY' + 1i*ZZ';
        end
    end
    tfr_noise_hard_nr(:,:,j) = tfr_noise_hard;
    tfr_noise_soft_nr(:,:,j) = tfr_noise_soft; 
    modes_hard(:,j) = 2*real(itfrstft_three_case_down(tfr_noise_hard,2,N,h,0));
    modes_soft(:,j) = 2*real(itfrstft_three_case_down(tfr_noise_soft,2,N,h,0));
end
end
