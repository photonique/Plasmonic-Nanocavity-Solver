%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
clear all; close all;
addpath("../")
aperiodic_depthsolver_scaled(2,[150 300],324,40,[30 40])
asymmetric_depthsolver_scaled(2,[300, 150, 150, 300],324,40,[ 40 30  30 40])
asymmetric_depthsolver_scaled_ttotal(2,[300, 150, 150, 300],324,40,[ 40 30  30 40])

## when the going gets tough
## the tough get going.
