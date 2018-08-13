%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%
% Copyright (C) 2009 Muthiah Annamalai
% 
% 02/23/09 : Wrapper for the intensity_far_field script, to modify it for
%            the aperiodic case, which allows it to be invoked similar to
%            all other aperiodic scripts in the genetic-alg project.
%
function Ifar_R_theta=aperiodic_intensity_far_field(Ex,theta_in,R,K,N,D,A,Ao)
  Ifar_R_theta = intensity_far_field(Ex,theta_in,R,K,N,[fliplr(D),0,D],A,Ao);
  return
end
