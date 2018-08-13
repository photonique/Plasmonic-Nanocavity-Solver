%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
## 
## 10/13/08: Apply  cloning and mutation to same population sets.
##           Remove Npop*2 in random population to Npop random elements.
##
## 10/10/08: Cloning and mutation applied to the depth parameter.
##           Remove the termination condition. We want the best!
## 
## 12/19/07: search in the aperiodic case, using the GA at 
##           a given number of corrugations only, since we cannot
##           keep track of the various size corrugations and various
##           corrugation periods (D) for each number of
##           corrugations ( N ).
## 
## 09/21/08: Add beam focus data to population tables. waist
##           is easily calculated as a0=1400*lambda/(560).
##           Add duplicate removal at each evolutionary stage.
## 
function population = ga_aperiodic_scaled(Ngen,xover_prob,optimize_choice,N,best_solution)

## FIXME: Add similarity removal in selection-culling process.

if ( nargin < 4 ), print_usage(), end;
if ( nargin != 5 )
  best_solution = [];
end
fprintf("GA: N = %g, Ngen = %g, xover_prob = %g, optimize = %s\n",N,Ngen,xover_prob,optimize_choice);
   
## optimization target.
optimize_T = true;
##optimize_T = false;
optimize_TL = ~ optimize_T;

Npop =100;    #keep this even #
mutate_prob = 0.1;
clone_rate = 1;
UNFIT_MAX = 1e3;

optimize_T = false;
optimize_TL = false;
optimize_global = false;

switch ( optimize_choice )
 case {'T'}
   disp('optimizing T')
   optimize_T = true;
 case {'TL'}
   disp('optimizing T/L')
   optimize_TL = true;
 case {'g'}
   disp('optimizing Global')
   optimize_global = true;
 otherwise
   error("Choice not understood");
end

fitness_T = @(X) ( 1/(1-X.T) );
fitness_TL = @(X) ( 1/(4-X.TL) ); ## just anything more than 1.13
fitness_global = @(X) ( 1/(1-X.T) + 1/(4-X.TL) );##emphasize T/L more than T.

if ( optimize_T ), fitness = fitness_T; fitness_max =  3000; end;
if ( optimize_TL ),  fitness = fitness_TL; fitness_max = 300; end;
if ( optimize_global), fitness = fitness_global; fitness_max = 3300; end

 ## fixed parameter a=40nm is implicit to solver.
 problem = @(N,D,Lambda,h) ( aperiodic_solver_scaled(N,D,Lambda, 40,h) );

 ## various alleles & their bit allocations
 col_T = 1;
 col_TL = 2;## 10bit
 col_lambda = 3; ## 10bit
 col_h = 4; ## 8bit
 col_fit = 5;
 col_prob = 6;
 col_focus = 7;
 col_d = 8; #col_d = 8:8+n-1 cols are filled
 
 col_total = col_d + N - 1;
 population = zeros(Npop,col_total);
 population(:,col_lambda) = randint( Npop, 1, 1023, 200);
 population(:,col_h) = randint( Npop, 1, 255, 40);
 population(:,col_fit) = zeros( Npop, 1);
 population(:,col_prob) = rand( Npop, 1);
 population(:,col_d:col_d+N-1) = randint( Npop, N , 200+N*40, 200);

 ## consolidate the best solutions 0-15, 16 is out of bounds.
 best_manual_solutions = zeros(15,col_total);
 best_manual_solutions(:,col_d:col_d+N-1) = randint(15,N,200+N*40,200);
 best_manual_solutions(:,col_lambda) = [280, 250, 240, 235, 230, 230, ...
		ones(1,9)*225];
 best_manual_solutions(:,col_h) = randint(15, 1, 255, 40);
 
 tl_best_sol = zeros(15,col_total);
 tl_best_sol(:,col_h) = randint( 15, 1, 255, 40);
 tl_best_sol(:,col_d:col_d+N-1) = repmat([500, 750, 800, 800, 850, 800, 750*ones(1,9)]',[1 N]);
 tl_best_sol(:,col_lambda) = [730, 960, 965, 945, 985, 920, 855, 855, 850*ones(1,7)];

 population = [best_solution; tl_best_sol; best_manual_solutions ; population];

 ## evaluate fitness.
 for idx = 1:rows(population)
      ##analp(population)
      val = population(idx,:);
      d = val(col_d:col_d+N-1); lambda = val(col_lambda); h = val(col_h);
      # eval fitness.
      yy = problem( N, d, lambda,h ); fit = fitness( yy );
      population( idx, col_h ) = yy.h;
      population( idx,col_fit) = fit;
      population( idx,col_T) = yy.T;
      population( idx,col_TL) = yy.TL;
      population( idx,col_focus) = yy.z0;
 end
 
 ## calc survival probabilities.
 ## probability of survival
 population(:,col_prob) =   population(:,col_fit)./sum( population(:,col_fit) );

 ## terminate least fit / culling 
 [val,idx] = sort(population(:,col_prob),'descend');
 population(idx(Npop+1:end),:)=[];

 unfit_population = []; #collects ancestors
 gen = 1;
 while gen < Ngen
	
	 ## re-calc survival probabilities.  
	 population(:,col_prob) = population(:,col_fit)./sum( population(:,col_fit) );

	## resort
	[val,idx] = sort(population(:,col_prob),'descend'); 
        fprintf('Fittest member @ gen %g\n',gen);
	val=population(1,:);
	val
        
	## clone top-1-fit of this population.
	cloned_pop = population(idx(1:clone_rate),:);

	## roulette wheel selection
	Nlarge=100*rows(population);
	prob_survival =  \
		population(:,col_fit)./sum( population(:,col_fit) );

	fitness_seq = weightedsel( prob_survival, Nlarge);
	
	## mutate population members at one location in pool,
	## by their p.d.f
	mutate_idx = randint(round(Npop*mutate_prob),1,Nlarge,1,true);
	mutate_idx = unique( fitness_seq( mutate_idx ) );

	mutated_pop = population(mutate_idx,:);
    for idx=0:N-1
	  mutated_pop( : , col_d + idx ) = ga_mutate_binvec( mutated_pop( :, col_d + idx ) );
    end
	mutated_pop( : , col_lambda ) = ga_mutate_binvec( mutated_pop( :, col_lambda ) );
	mutated_pop( : , col_h ) = ga_mutate_binvec( mutated_pop( :, col_h ) );

	xover_pop = mutated_pop;
    for idx=0:N-1
    	xover_pop( : , col_d + idx ) = ga_xover_binvec( xover_pop( :, col_d + idx ) );
    end
	xover_pop( : , col_lambda ) = ga_xover_binvec( xover_pop( :, col_lambda ) );
	xover_pop( : , col_h ) = ga_xover_binvec( xover_pop( :, col_h ) );
    
	## evaluate fitness of xover+mutated and cloned population
	new_pop = [ xover_pop; cloned_pop];

	## remove n =0 case, a known minimum
	## remove illegal data point solutions
	## remove (lambda, d) < 200 as well as h < 30.
	idx = find( new_pop(:,col_lambda) < 200 );	
	new_pop(idx,:) = [];
	idx = find( new_pop(:,col_h) < 30 );
	new_pop(idx,:) = [];

	for idx = col_d:col_d+N-1
            idy = find( new_pop(:,idx) < 100 );
	    new_pop(idy,:) = [];
	end


    for idx = 1:size(new_pop,1)
          val = new_pop(idx,:);
          d = val(col_d:col_d+N-1); lambda = val(col_lambda); h = val(col_h);
          # eval fitness.
          yy = problem( N, d, lambda,h ); fit = fitness( yy );
          new_pop( idx, col_h ) = yy.h;
          new_pop( idx,col_fit) = fit;
          new_pop( idx,col_T) = yy.T;
          new_pop( idx,col_TL) = yy.TL;
     end
    
    ## consolidate population
	population = [population; new_pop ];
    
	## calc survival probabilities.
	tot = sum( population(:,col_fit) );
    prob_survival = population(:,col_fit)./tot;
    population(:,col_prob) = prob_survival;

    ## terminate least fit / culling 
    [val,idx] = sort(population(:,col_prob),'descend');
    unfit_population = [unfit_population; population(idx(Npop+1:end),:)];
    population(idx(Npop+1:end),:)=[];

    ## keep only the top UNFIT_MAX unfit and in descending order.
    if ( rows( unfit_population ) > UNFIT_MAX )
      [val,idx] = sort(unfit_population(:,col_prob),'descend');
      unfit_population = unfit_population(idx(1:UNFIT_MAX),:);
    end
    
    ## remove duplicates in the population and fill with topmost 
    ## unused populations.
    [upop,idpop,idupop]=unique(population,'rows');
    population = upop;
    idxmissing = Npop - rows(population);
    if ( idxmissing < rows(unfit_population) )
        population = [population; unfit_population(1:idxmissing,:)];
        unfit_population=unfit_population(idxmissing+1:end,:);
    end
    
    [val,idx] = sort(population(:,col_prob),'descend');
    population = population(idx,:);
    
    gen = gen + 1;
 end

 ## re-calc survival probabilities.
 prob_survival = population(:,col_fit)./sum( population(:,col_fit) );
 population(:,col_prob) = prob_survival;

 ## resort and retain fittest.
 [val,idx] = sort(population(:,col_prob),'descend');

 ## topmost is the best now
 population = population(idx,:);
 
 ## display
 fprintf('generation = %d, best-solution\n',gen);
 fprintf('%g \n',population(1,:));
end
