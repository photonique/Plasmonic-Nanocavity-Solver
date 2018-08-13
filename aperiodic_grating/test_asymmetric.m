%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## 07/16/09 : Sanity checks for the Asymmetric solver
## 
clear all; close all; addpath("../")
n = 3; D = [ 133, 313, 304 ]; h = 61; lambda = 336; a = 40;
Dasym = [ fliplr( D), D ];
fprintf(" expected T = 0.318, T/L = 0.57 \n");

d = mean( D ); 
#D = repmat( d, [ 1, 3]);  Dasym = repmat( d, [ 1, 6 ] );
fprintf(" ******* PERIODIC & SYMMETRIC ********* \n")
T_TL_max=solver_scaled_nanocavity( n, d, lambda, a )

fprintf(" ******* APERIODIC & SYMMETRIC ********* \n")
T_TL_max = aperiodic_solver_scaled( n, D, lambda, a, h )

#Dasym = [ fliplr(D) - 2, D ];
#Dasym = [ D, D ];
fprintf(" ******* APERIODIC & ASYMMETRIC ********* \n")
T_TL_max = asymmetric_solver_scaled( n, Dasym, lambda, a, h )

