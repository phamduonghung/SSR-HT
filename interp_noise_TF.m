function [tfr_interp] = interp_noise_TF(cas,tfr_noise,Band_noise)
 %extrapolation of the noise from the data outside the region of interest 
 B = size(tfr_noise);

 %the coordinates
 X = ones(B(1),1)*(1:B(2));
 Y = ((0:B(1)-1)/B(1))'*ones(1,B(2));
 
 %the region ones keeps
 A =zeros(B);
 A(:,:) = Band_noise;
 if cas == 1
  %2D interpolation  
  %the corresponding X and Y, columnwise
  XA = X(A == 1);
  YA = Y(A == 1);
  ZA = tfr_noise(A==1);
  Real_ZA = real(ZA);
  Imag_ZA = imag(ZA);
  
  [x1,y1] = ndgrid((0:B(1)-1)/B(1),1:B(2));
 
  z1 = griddata(XA,YA,Real_ZA,y1(:),x1(:),'cubic');
  z2 = griddata(XA,YA,Imag_ZA,y1(:),x1(:),'cubic');
  tfr_interp = reshape(z1+1i*z2,B(1),B(2));
 else
  %one-dimensional interpolation of the real part
  tfr_interp = zeros(B);
  Z = real(tfr_noise).*(A==1);
  ZZ = imag(tfr_noise).*(A==1);
  C = Y.*(A==1);
  for k=1:B(2),
   Val = Z(:,k);
   pos = C(:,k);
   Val1 = Val(Val ~= 0);
   pos  = pos(Val ~= 0);
   yy  = spline(pos,Val1,(0:B(1)-1)/B(1));
   Val = ZZ(:,k);
   pos = C(:,k);
   Val1 = Val(Val ~= 0);
   pos  = pos(Val ~= 0);
   yy1  = spline(pos,Val1,(0:B(1)-1)/B(1));
   tfr_interp(:,k) = yy+1i*yy1;
  end
 end 
end
