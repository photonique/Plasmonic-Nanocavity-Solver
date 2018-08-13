%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C) 2009 Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
##
## 02/03/09 : add a priming data set option.
## 
## 01/30/09 : 
## 
## Generic Particle Swarm Optimization (PSO) for a maximization problem,
## given variables, their ranges and the PSO parameters.
## 
## problem: function handle takes Ndims arguments, nothing optional or extra,
## and returns a single number proportional to the fitness. problem must
## take the argument as a [ 1 x Ndims ] col vector. Optionally can take an
## extras argument if supplied.
## 
## Ndims : size of the decision / solution space (i.e R^3, R^Ndims etc),
## where R is Riemann set.
## 
## RINFO: is [ Ndims x 2 ] matrix with 1-col of start pt, 2-col of end pt.
## range of the variables: eg: r_start = -100; r_end = +100;, for 1D problem.
## 
## FINFO: fitness information
## min_fit = 0, max_fit = 1;
## 
## PSO_PARAMS: control variables for PSO [ 1 x 7 ] col-vector.
## -inertial weights of PSO, cognitive , social accel factor, max vel.
## WF = 1; C1 = 2; C2 = 2; max_velocity = 10;
## -generations, population size, tolerance in % from max-fit.
## Ngen = 200; Npop = 150; Tol = 1;
## 
## Allow extras, arguments to pass on to the problem. It is superfluous,
## but makes things easier.
##
## Add noeval_flag, to not evaluate the population member that has
## positions outside allowed hyperspace. Correct implementation of 
## checking range limits in dimensions of hyperspace.
## 
## Use MAX_VELOCITY parameter, as not scaled with fitness.
## 
function [ population, gbest_pos, gbest_fitness, gen ] = pso( problem, Ndims, RINFO, FINFO, PSO_PARAMS
							     ,EXTRAS , best_sol )

if ( nargin < 4 || nargin > 7 )
  print_usage ( ) ;
end

## do some sanity check on dimensions of the arguments.
if ( size( RINFO ) != [ Ndims, 2 ] )
  error(' RINFO : need to be of size [ Ndims, 2 ]')
elseif ( size( FINFO ) != [ 1, 2 ] )
  error(' FINFO : need to be of size [ 1, 2 ]')
elseif ( nargin >= 5 )
  if ( size( PSO_PARAMS ) != [ 1, 7 ] )
	error(' PSO_PARAMS : need to be of size [ 1, 7 ]')
  end
end

## unpack the parameters & setup the problem.
r_start = RINFO( :, 1 )'; r_end = RINFO ( :, 2 )';
min_fit = FINFO( 1 ); max_fit = FINFO ( 2 );

r_range = (r_end - r_start);

## makes a coordinate in N-space between [ r_start to r_end ].
make_coord = @( X ) ( r_start + r_range.*X );

## PSO parameters
if ( nargin >= 5 )
  WF = PSO_PARAMS( 1 ); C1 = PSO_PARAMS( 2 ); 
  C2 = PSO_PARAMS( 3 ); MAX_VELOCITY = PSO_PARAMS( 4 );
  Ngen = PSO_PARAMS( 5 ); Npop = PSO_PARAMS( 6 );
  Tol = PSO_PARAMS( 7 );
else
  ## ok defaults.
  WF = 1; C1 = 2; C2 = 2; MAX_VELOCITY = 10;
  Ngen = 200; Npop = 150; Tol = 1;
end

if ( nargin < 6 )
  EXTRAS = [];
end

if ( nargin < 7 )
  best_sol = [];    
end

population = { };
fitness = [];

## initial assignment of position, velocity,
## and fitness evaluation round
gbest_fitness = min_fit; gbest_pos = [];
Npop = rows( best_sol ) + Npop ;
for idx = 1 : Npop

  if ( idx <= rows( best_sol ) )
    X = best_sol( idx, : );
  else    
    X =  make_coord( rand(1,Ndims) );
  end
  
  ## current values are implied here,
  population{ idx }.pos = X;
  population{ idx }.fitness = problem( X, EXTRAS );

  population{ idx }.pbest_pos = population{ idx }.pos;
  population{ idx }.pbest_fitness = population{ idx }.fitness;

  population{ idx }.velocity = (2*rand(1,Ndims) - 1)*MAX_VELOCITY;

  if ( gbest_fitness < population{ idx }.fitness )
    gbest_pos = population{ idx }.pos;
    gbest_fitness = population{ idx }.fitness;
  end

end


## begin the PSO loop
gen = 1;
while ( gbest_fitness < max_fit*(1-Tol/100) && gen < Ngen )
  
  gen = gen + 1;
  
  ## fitness proportional max velocity, generally slow down as you reach
  ## top, and do finer/granular searches.
  ## max_velocity = MAX_VELOCITY*( max_fit - gbest_fitness )/max_fit ;

  fprintf('At gen %g best population member : fitness = [ %g]\n', gen, gbest_fitness);

  ## record
  fitness( gen ) = gbest_fitness;
  
  ## 
  ## calculate the particle velocity and position
  ## and update the fitness.
  ## 
  for idx = 1 : Npop

    ## update velocity
    population{ idx }.velocity = WF*population{ idx }.velocity + \
	C1*( 2*rand(1,Ndims) - 1 ).*(population{ idx }.pbest_pos - population{ idx }.pos) + \
	C2*( 2*rand(1,Ndims) - 1 ).*( gbest_pos  - population{ idx }.pos) ;

    noeval_flag = zeros( Npop, 1 );

    ## clamp coordinates in each dimension; when exceed range use
    ## velocity change ( reflecting wall ) concepts.
    for dim_idx = 1 : Ndims

     if ( abs( population{ idx }.velocity( dim_idx ) ) >= MAX_VELOCITY )
       population{ idx }.velocity( dim_idx ) = -sign( population{ idx }.velocity( dim_idx ) )*MAX_VELOCITY ;
       noeval_flag( idx ) = true;
     end

     pos = population{ idx }.pos;

     if ( pos( dim_idx ) >= r_end( dim_idx ) || pos( dim_idx ) <= r_start( dim_idx ) )
       population{ idx }.velocity( dim_idx ) = -population{ idx  }.velocity(  dim_idx );
       noeval_flag( idx ) = true;
     end

    end
    
    
    ## update pos
    X = population{ idx }.pos;
    X = X + population{ idx }.velocity;

    population{ idx }.pos = X;

    ## evaluate 
    if ( ~ noeval_flag( idx ) )
	  X = population{ idx }.pos;
	  population{ idx }.fitness =  problem( X, EXTRAS );
    end

    ## re-assign
    population{ idx }.pbest_fitness = max( population{ idx }.pbest_fitness, population{ idx }.fitness );
    
    if ( gbest_fitness < population{ idx }.fitness )
      gbest_fitness = population{ idx }.fitness;
      gbest_pos = population{ idx }.pos;
    end

  end
  
end

fprintf('At gen %g best population member : fitness = [ %g]\n', gen, gbest_fitness);

return;
end
