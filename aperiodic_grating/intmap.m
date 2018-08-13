%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 08/12/09
## Check effect of the aperiodic sequences on the maps.
## 
## 08/11/09
## Asymmetric and aperiodic nanocavity with two arms
## optimized for the left and right tuned to periodic
## case N=2, resonant with the wavelengths  238 and 730.
## 
close all; clear all; more off;
addpath("../")
do_intensity_map = true;
asymmetric = true;

N = 5;
a = 40;
## D = repmat(d,[1, N]);

if ( ! asymmetric )
  D = [ 290 290];
  H = 90.1;
else
  D = [ 290 290 718 718 ];
  H = [ 90.1 90.1 156.1 156.1 ];

  ## trim thue-morse to correct length
  apseq = thue_morse1d(N+(N<=3)*2,1);
  apseq = fibo1d(2*N,1);
  
  L = round(length(apseq)/2)*2;
  apseq = [apseq(L/2-floor(N)+1:L/2+ceil(N))];

  #include symmetry ??
  #apseq = [fliplr(apseq), apseq];

  D = [290, 718](apseq);
  H = [90.1, 156.1](apseq);
  
  D
  H
end
d = mean(D);
h = mean(H);

theta_in = -90:+90; R = 1; Ao = 1;

wavelengths = 200:10:1000;

Npts = length( wavelengths );
T = zeros(1,Npts); TL = zeros(1,Npts); Ttotal = zeros(1,Npts);

for idx = 1:Npts
  
  Lambda = wavelengths( idx );
  
  if ( do_intensity_map )
    if ( asymmetric )
      [Ttotal( idx ),Ex] = ttotal_asymmetric_solver(N,D,Lambda,a,H);
      K = 2*pi/Lambda;
      Ifar(idx,:) = asymmetric_intensity_far_field(Ex,theta_in,R,K,N,D,a,Ao);
    else
      [Ttotal( idx ),Ex] = ttotal_periodic_solver(N,d,Lambda,a,H);
      K = 2*pi/Lambda;
      Ifar(idx,:) = aperiodic_intensity_far_field(Ex,theta_in,R,K,N,D,a,Ao);    
    end
  else
    yy = aperiodic_solver_scaled(N,D,Lambda, a,H );
    T( idx ) = yy.T;
    TL( idx ) = yy.TL;
    Ttotal( idx ) = yy.T*(1 + 1/yy.TL);
  end
  
end

geom_info_str = sprintf(" N = %d, H = %g, D = %g, A = %g",N,h,d,a);

if ( do_intensity_map )

  mesh(theta_in,wavelengths,Ifar);
  set(gcf(),'color','white')
  xlabel('Azimuthal Angle \theta (degree)')
  axis([min(theta_in)-1,1+ max(theta_in)])
  ylabel('\lambda (nm)')
  shading('interp')
  colorbar()
  view(0,90)
  title(['Intensity(\lambda,\theta),',geom_info_str]);

else

  plot( wavelengths, T*100, '-ob')
  set(gcf(),'color','white')
  title(['Transmission spectra : Coupling to Gaussian mode',geom_info_str])
  xlabel('Wavelength [nm]')
  ylabel("Coupling to Gaussian mode [%]")

  figure()
  plot( wavelengths, TL, '-og')
  set(gcf(),'color','white')
  title(['Transmission spectra : T/L ',geom_info_str])
  xlabel('Wavelength [nm]')
  ylabel(" T/L ratio")

  figure()
  plot( wavelengths, Ttotal*100, '-or')
  set(gcf(),'color','white')
  title(['Total transmission spectra',geom_info_str])
  xlabel('Wavelength [nm]')
  ylabel("Total transmission [%]")

endif
