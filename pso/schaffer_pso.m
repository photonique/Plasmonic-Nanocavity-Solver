%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## 01/29/09 : 
## 

clear all
close all
more off

do_anim = false;

fitness = @( fit ) 1/( fit );
problem = @( x, y )  fitness( schaffer_f6( [x, y] ) );

## range of the variables
r_start = -100; r_end = +200; r_range = (r_end - r_start);
min_fit = -1; max_fit = 1;

## X is in size 0 - 1.
make_coord = @( X ) ( r_start + r_range.*X );

## PSO parameters
WF = 1; C1 = 2; C2 = 2; max_velocity = 0.15;

## size of population, no of generations
Npop = 100; Ngen = 50;

population = { };
fitness = [];
## 
## each population member ( bee ) has the following fields
## fitness
## pos = [ x, y]; and pbest_pos.
## velocity
## similarly one global copy of  gbest_pos, gbest_fitness
## 

## initial assignment of position, velocity,
## and fitness evaluation round
gbest_fitness = min_fit; gbest_pos = [];
for idx = 1 : Npop
  X =  make_coord( rand(1,2) );

  ## current values are implied here,
  population{ idx }.pos = [ X(1), X(2) ]; ##emphasize 2-elem nature
  population{ idx }.fitness = problem( X(1), X(2) );

  population{ idx }.pbest_pos = population{ idx }.pos;
  population{ idx }.pbest_fitness = population{ idx }.fitness;

  population{ idx }.velocity = (2*rand(1,2)-1)*max_velocity;

  if ( gbest_fitness < population{ idx }.fitness )
    gbest_pos = population{ idx }.pos;
    gbest_fitness = population{ idx }.fitness;
  end

end


## begin the PSO loop
for gen = 1 : Ngen

  fprintf('At gen %g best population member : coord [ %g, %g ], fitness = [ %g]\n',\
	                                  gen,
					  gbest_pos(1),
					  gbest_pos(2),
					  gbest_fitness);

  ## record
  fitness( gen ) = gbest_fitness;
  
  ## 
  ## calculate the particle velocity and position
  ## and update the fitness.
  ## 
  POS = zeros( Npop ,2 );
  for idx = 1 : Npop

    ## update velocity
    population{ idx }.velocity = WF*population{ idx }.velocity + \
	C1*(2*rand(1,2)-1).*(population{ idx }.pbest_pos - population{ idx }.pos) + \
	C2*(2*rand(1,2)-1).*( gbest_pos  - population{ idx }.pos) ;

    ## clamp velocity using reflecting walls
    if ( abs( population{ idx }.velocity ) > max_velocity )
      population{ idx }.velocity = -sign( population{ idx }.velocity )*max_velocity ;
    end

    ## update pos
    X = population{ idx }.pos;
    X = X + population{ idx }.velocity;

    population{ idx }.pos = X;

    ## evaluate 
    X = population{ idx }.pos;
    POS( idx, : ) = X;
    population{ idx }.fitness =  problem( X(1), X(2) );
    
    ## re-assign
    population{ idx }.pbest_fitness = max( population{ idx }.pbest_fitness, population{ idx }.fitness );
    
    if ( gbest_fitness < population{ idx }.fitness )
      gbest_fitness = population{ idx }.fitness;
      gbest_pos = population{ idx }.pos;
    end

  end
  
  ## display the swarm
  plot( POS( :, 1 ), POS( :, 2), sprintf('ob; swarm @ gen %g;',gen) );
  title(sprintf('Swarm Population @ gen %g;',gen));
  axis( [ -2 2 -2 2 ] );
  xlabel(' solution X2'); ylabel(' solution X1');
if ( do_anim )
  print('-dpng',sprintf('swarm_opt_%d.png',gen))
end

end

fprintf('At gen %g best population member : coord [ %g, %g ], fitness = [ %g]\n',\
	                                  Ngen,
					  gbest_pos(1),
					  gbest_pos(2),
					  gbest_fitness);

plot(1:Ngen,1./fitness,'-o;best fitness as function of generations of PSO;')
title( ' PSO fitness grow with generations ')


if ( do_anim )
 system("convert -delay 50 -loop 1 swarm_opt_*.png PSO.gif")
 system("animate PSO.gif")
end
