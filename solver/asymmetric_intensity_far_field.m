%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%
% Copyright (C) 2009 Muthiah Annamalai
% 
% 07/23/09 : Wrapper for the intensity_far_field script, to modify it for
%            the asymmetric case, for argument order similar to other scripts.
%            Useful to visualize the I(\lambda,\theta) plots in 2D.
%
function Ifar_R_theta=asymmetric_intensity_far_field(Ex,theta_in,R,K,N,D,A,Ao)
  Ifar_R_theta = intensity_far_field(Ex,theta_in,R,K,N,[D(1:N),0,D(N+1:2*N)],A,Ao);
  return
end
