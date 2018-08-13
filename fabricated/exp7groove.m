%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## 02/23/09 : model for the experimental piece fabricated by Lei Zhu.
## Periodic nanocavity
## 
close all; clear all; more off;
addpath("../")
do_intensity_map = true;
asymmetric = !true;

N = 7;
a = 70;
d = 500;
h = 110;
## D = repmat(d,[1, N]);

if ( ! asymmetric )
  D = d*ones(1,N);
else
  D = d + ( 20*(2*(rand(1,2*N)>0.5) - 1) );
end


theta_in = -90:+90; R = 1; Ao = 1;

wavelengths = 200:15:1000;

Npts = length( wavelengths );
T = zeros(1,Npts); TL = zeros(1,Npts); Ttotal = zeros(1,Npts);

for idx = 1:Npts
  
  Lambda = wavelengths( idx );
  
  if ( do_intensity_map )
    if ( asymmetric )
      [Ttotal( idx ),Ex] = ttotal_asymmetric_solver(N,D,Lambda,a,h);
      K = 2*pi/Lambda;
      Ifar(idx,:) = asymmetric_intensity_far_field(Ex,theta_in,R,K,N,D,a,Ao);
    else
      [Ttotal( idx ),Ex] = ttotal_periodic_solver(N,d,Lambda,a,h);
      K = 2*pi/Lambda;
      Ifar(idx,:) = aperiodic_intensity_far_field(Ex,theta_in,R,K,N,D,a,Ao);    
    end
  else
    yy = aperiodic_solver_scaled(N,D,Lambda, a,h );
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
