%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 01/30/09 : Solve the same sinc2d problem using PSO, albeit the
##            generic version, for testing the PSO script
clear all
close all
more off

Ndims = 2;
problem = @( X ) sinc( X(1)/pi ).*sinc( X(2)/pi );
RINFO = [ -10 +10;
	  -10 +10 ];
FINFO = [ 0 1 ];
WF = 1; C1 = 2; C2 = 2; MAX_VELOCITY = 3;
Ngen = 200; Npop = 150;50; Tol = 1;
PSO_PARAMS = [ WF, C1, C2, MAX_VELOCITY, Ngen, Npop, Tol ];
rand( "seed",  time( ) )
[ population, gbest_pos, gbest_fitness, gen ] = pso( problem, Ndims, RINFO, FINFO, PSO_PARAMS )
