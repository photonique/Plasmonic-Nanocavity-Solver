%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
##
## 01/30/09 : Broaden the search space. Add restriction operator in
##            dimensions, and velocity clamping, also based on
##            coordinate ranges. 2*rand(x,y)-1 choice of random number is
##            very highly convergent to the best solution ~ 0.98 for the
##            cases of our large search space -100 to +100, over 200 generations.
##
##            With an increased population 10x more than earlier (150, w.r.t 15) 
##            you can clearly see the surety in the search space
##            convergen in about 10 generations! Amazing. All this
##            without fancy inertial weights, but with max velocity tapering
##            proportional to fitness. Also most "finds" occur <= 400 generations!!
## 
## 01/28/09 : Find the maxima of the 2D SINC function (0,0) point in the
##            decision space, for the objective function value of 1, using
##            the Particle Swarm Optimization ( PSO ).
##            See: "Particle Swarm Optimization", Kennedy, Eberhart 1995.
##            
##            Wow! PSO rocks. Resolution of search space is governed by
##            choice of maximum velocity. Exercise choice while picking this.
##            
##            Eliminate the need to broadcast the global best values,
##            and keep Npop-instances of the global pos, fitness.
## 
##            PSO global fitness and position need not be always present
##            in the population. It is just a record, and populations generally
##            overshoot the global maxima in a typical problem.
##            
##            Make random numbers in the velocity calculation go from -1:+1.
##            This accelerates convergence & shows the visual "swarming". 
##            Add default false animation flag.

clear all
close all
more off

do_anim = false;
rand("seed",time());

problem = @( x, y) sinc(x/pi).*sinc(y/pi) ;

my_rand = @( x, y ) 2*rand(x,y)-1 ;


## range of the variables
r_start = -100; r_end = +100; r_range = (r_end - r_start);
min_fit = 0; max_fit = 1;

## X is in size r_start to r_end
make_coord = @( X ) ( r_start + r_range.*X );

## PSO parameters
WF = 1; C1 = 2; C2 = 2; max_velocity = 10;

## size of population, no of generations
Npop = 150; Ngen = 400;

## choice of velocity ranges
## velocity_range = linspace(10,0.01,Ngen);
## choice of inertial weights.
## WF = linspace(0.9,0.1,Ngen);


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

  population{ idx }.velocity = (my_rand(1,2))*max_velocity;

  if ( gbest_fitness < population{ idx }.fitness )
    gbest_pos = population{ idx }.pos;
    gbest_fitness = population{ idx }.fitness;
  end

end


## begin the PSO loop
##while ( gbest_fitness < 0.99 ) 
##   gen = gen + 1;
## for gen = 1 : Ngen
gen = 1;
while ( gbest_fitness < 0.99  && gen < Ngen )
  gen = gen + 1;
  
  ## fitness proportional max velocity, generally slow down as you reach
  ## top, and do finer/granular searches.
  max_velocity = 10*( 1 - gbest_fitness ) ;

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
	C1*( my_rand(1,2) ).*(population{ idx }.pbest_pos - population{ idx }.pos) + \
	C2*( my_rand(1,2) ).*( gbest_pos  - population{ idx }.pos) ;

    ## clamp coordinates in each dimension; when exceed range use
    ## velocity change ( reflecting wall ) concepts.
    if ( abs( population{ idx }.velocity(1) ) >= max_velocity )
      population{ idx }.velocity(1) = -sign( population{ idx }.velocity(1) )*max_velocity ;
    end
    if ( abs( population{ idx }.velocity(2) ) >= max_velocity )
      population{ idx }.velocity(2) = -sign( population{ idx }.velocity(2) )*max_velocity ;
    end
    pos = population{ idx }.pos;
    if ( pos(1) >= r_end || pos(1) <= r_start )
      population{ idx }.velocity(1) = -population{ idx }.velocity(1);
    end    
    if ( pos(2) >= r_end || pos(2) <= r_start )
      population{ idx }.velocity(2) = -population{ idx }.velocity(2);
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
  if ( do_anim )
    plot( POS( :, 1 ), POS( :, 2), sprintf('ob; swarm @ gen %g;',gen) );
    title(sprintf('Swarm Population @ gen %g;',gen));
    axis( [ r_start, r_end, r_start, r_end  ] );
    xlabel(' solution X2'); ylabel(' solution X1');
    print('-dpng',sprintf('swarm_opt_%d.png',gen));
  end

end

fprintf('At gen %g best population member : coord [ %g, %g ], fitness = [ %g]\n',\
	                                  gen,
					  gbest_pos(1),
					  gbest_pos(2),
					  gbest_fitness);

if ( do_anim )
 plot(1:Ngen,fitness,'-o;best fitness as function of generations of PSO;')
 title( ' PSO fitness grow with generations ')
 system("convert -delay 50 -loop 1 swarm_opt_*.png PSO.gif")
 system("animate PSO.gif")
end
