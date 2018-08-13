%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%% generate figures for the JMO article 03/26/10
%% 
%% All designs for lambda =560nm resonance.
%% 
clear all
close all
%% Aperiodic spacing optimization with uniform groove depth
N=6; D = [218,503,510,534,162,324]; H = 101*ones(6,1); A =40;
figure
drawcavity(N,D,H,A)

%% same with nonuniform groove depth
N=8; D = [239,120,341,517,516,533,558,130]; 
     H = [103,97,102,95,104,105,106,119]; A = 40;
figure
drawcavity(N,D,H,A)

%% the ultimate best cavity design at Lambda = 560
N=9; D = [238 , 145 , 331 , 152 , 387 , 455 , 60 , 514 , 508]; 
     H = [99 , 95 , 100 , 94 , 100 , 12 , 100 , 105 , 112]; A = 40;
figure
drawcavity(N,D,H,A)
