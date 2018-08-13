%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%%
## 
## 08/6/09 : Sanity checks for the Asymmetric solver
## 
clear all; close all; addpath("../")
n = 5;  h = 61; lambda = 336; a = 40;

for iters = 1:1000
 D = [ 200, 400, 133, 313, 304 ] + 20*(2*(rand(1,n)>0.5) - 1);
 Dasym = [ fliplr( D), D ];
 
 ## fprintf(" ******* APERIODIC & SYMMETRIC ********* \n")
 T_TL_aper = aperiodic_solver_scaled( n, D, lambda, a, h );

 ## fprintf(" ******* APERIODIC & ASYMMETRIC ********* \n")
 T_TL_asym = asymmetric_solver_scaled( n, Dasym, lambda, a, h );

 fprintf(" [ AP| %g , ASYM| %g ]\n",T_TL_aper.T, T_TL_asym.T);
 
 assert( abs(T_TL_aper.T - T_TL_asym.T ) < 0.01 )
 assert( abs(T_TL_aper.TL - T_TL_asym.TL ) < 0.01 )
 
end
