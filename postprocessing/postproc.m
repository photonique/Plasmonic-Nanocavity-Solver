%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## 02/27/09 : calculate dispersion relations
## 
clear all 
more off

load "fab7groove-Ithetalambda.mat";
## theta_in = theta_in( find( abs(theta_in) < 20 ) );
## calculate 2 branches of lambda SPP
lambda_spp = zeros(rows(theta_in),2);
for i = 1:length(theta_in)
  if ( 0 )
    plot(wavelengths,Ifar(:,i),'-ob')
    title(["Spectrum at: angle ",num2str(theta_in(i))])
    xlabel("Wavelength")
    ylabel("Intensity")
  end
  
  [tmp,idx]=sort(Ifar(:,i),'ascend');
  fprintf("Minimum wavelength : = %g at Angle %g\n",wavelengths(idx(1)),theta_in(i))
  lambda_spp( i, 1:2 ) = wavelengths( [idx(1),idx(end)] );
  ##pause()
end

fprintf(" Theta : Lambda SPP 1, 2 \n")
[reshape(theta_in,[length(theta_in), 1]), lambda_spp]

k = 2*pi./(lambda_spp*1e-9)*1e-6; #k units of inverse-microns
G = (2*pi/500e-9)*1e-6;

## phase matching condition to launch a plasmon wave
## this is proper way to calculate the Kspp
kspp = G + k.*repmat(sin(reshape(theta_in,[length(theta_in),1])*pi/180),[1,2]);

w = 3e8./(lambda_spp*1e-9);
plot( kspp(:,1), w(:,1),'or;high energy branch;', kspp(:,2), w(:,2),'+g;low energy branch;')

#plot( 2*pi./(500*1e-9).*sin((theta_in)*pi/180), -w(:,1) , '-or' );
#plot( 2*pi./(lambda_spp(:,1)*1e-9).*sin((theta_in)*pi/180), -w(:,1) , '-ob' );
#plot( 2*pi./(500*1e-9).*sin((theta_in)*pi/180), -w(:,1) , '-or' );
#plot( kspp, w(:,1), '-ob')
