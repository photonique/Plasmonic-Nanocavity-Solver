%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C)  2007, 2008 Muthiah Annamalai <muthiah.annamalai@uta.edu>
##
##
## 10/29/08: Following advice from advisor, I have used a random +/- 1
##           unit steps for the cloning & mutation of the h-parameter in
##           GA. Earlier I used simple-GA also on h-parameter which
##           doesnt work. Clearly moving to random search or adjustment
##           weights works better.
##           
## 10/22/08: Donot evaluate best-solution first, assume its calculated.
##           Experiment with new fitness functions. Correct way to
##           eliminate illegal mutants in the population. Population crowds
##           way too easily! Remove elitism. Iteration to achieve proper 
##           legal mutation populations sizes . Convergence is really
##           slow once I removed elitism!
## 
## 10/18/08: Certain choice of parameters lend to singular matrices
##               and we necessarily need to change the solver for
##                inverting the matrix. Tag unsolvable cases where Gab is 
##		    badly conditioned as -1, and remove from population before
##                calculating the fitness. This preserves concept of GA moving
##                towards a fit population.
## 
## 10/10/08: Vary depth, wavelength and period of corrugations as well
##           in aperiodic-symmetric manner. Cloning and mutation applied to the depth parameter.
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
function population = ga_depth_aperiodic_scaled(Ngen,xover_prob,optimize_choice,N,best_solution)

## FIXME: Add similarity removal in selection-culling process.

if ( nargin < 4 ), print_usage(), end;
if ( nargin != 5 )
  best_solution = [];
end
fprintf("GA: N = %g, Ngen = %g, xover_prob = %g, optimize = %s\n",N,Ngen,xover_prob,optimize_choice);

## optimization target.
optimize_T = true;
optimize_TL = ~ optimize_T;

Npop = 50;    #keep this even #
mutate_prob = 0.1;
clone_rate = 1;
UNFIT_MAX = 1e2;

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

## fitness_T = @(X) ( 1000*X.T );
## fitness_TL = @(X) ( 100*X.TL ); ## just anything more than 1.13
## fitness_global = @(X) ( 100*X.T + 1000*X.TL );##emphasize T/L more than T.

##fitness_T = @(X) ( 1/(1-X.T) );
##fitness_TL = @(X) ( 1/(4-X.TL) ); ## just anything more than 1.13
##fitness_global = @(X) ( 1/(1-X.T) + 1/(4-X.TL) );##emphasize T/L more than T.

if ( 0 )
    fitness_T = @(X) ( 1 - (0.5-X.T)^2 );
    fitness_TL = @(X) ( 4^2 - (4-X.TL)^2 ); ## just anything more than
    fitness_global = @(X) ( 1 - (0.6-X.T)^2 + 4^2 - (4-X.TL).^2 );##emphasize T/L more than T.
end

if ( 0 )
    fitness_T = @(X) ( 1 - exp(-4*X.T.^2) );
    fitness_TL = @(X) ( 1 - exp(-4*X.TL.^2) ); ## just anything more than 1.13
    fitness_global = @(X) ( 1 - exp(-4*X.T) + 1 - exp(-4*X.TL.^2) );##emphasize T/L more than T.
end

 fitness_T = @(X) ( 1000*X.T );
 fitness_TL = @(X) ( 1000*X.TL ); ## just anything more than 1.13
 fitness_global = @(X) ( 1000*X.T + 1000*X.TL );##emphasize T/L more than T.


if ( optimize_T ), fitness = fitness_T; fitness_max =  Inf; end;
if ( optimize_TL ),  fitness = fitness_TL; fitness_max = Inf; end;
if ( optimize_global), fitness = fitness_global; fitness_max = Inf; end

## fixed parameter a=40nm is implicit to solver.
problem = @(N,D,Lambda,H) ( aperiodic_depthsolver_scaled(N,D,Lambda, 40,H) );

 ## various alleles & their bit allocations
 col_T = 1;
 col_TL = 2;## 10bit
 col_lambda = 3; ## 10bit
 col_focus = 4; ## 8bit
 col_fit = 5;
 col_prob = 6;
 col_d = 7; #col_d = 7:7+n-1 cols are filled
 col_h = 7+N; #col_d = 7+n:7+n+n-1 cols are filled
 
 col_total = col_h + N - 1;
 population = zeros(Npop,col_total);
 population(:,col_lambda) = randint( Npop, 1, 1023, 200);
 population(:,col_fit) = zeros( Npop, 1);
 population(:,col_prob) = rand( Npop, 1);
 population(:,col_d:col_d+N-1) = randint( Npop, N , 300+N*100, 200);
 %scaled to the wavelength of interest, unlike earlier independent of wavelength
 %lambda/4 has some singularity at this point in the calculations
 population(:,col_h:col_h+N-1) = round(diag(population(:,col_lambda)/4.1)*rand( Npop, N));
 %population(:,col_h:col_h+N-1) = randint( Npop, N,255,10);
 
 %% cleanup illegal population
    idx = find( population(:,col_lambda) < 200 );
    population(idx,:) = []; 

    rmlist = zeros(rows(population),1);
    for idx = col_d:col_d+N-1
      idy = find( population(:,idx) < 10 );
      rmlist(idy) = 1;
    end
    population(find(rmlist>0),:) = [];
    
    rmlist = zeros(rows(population),1);
    for idx = col_h:col_h+N-1
      idy = find( population(:,idx) < 10 );
      rmlist(idy) = 1;
    end
    population(find(rmlist>0),:) = [];

 %club with best sol assuming its already evaluated.
 population = [best_solution; population];

 ## evaluate fitness.
 for idx = 1:rows(population)
      ## analp(population)
      val = population(idx,:);
      d = val(col_d:col_d+N-1); lambda = val(col_lambda); h = val(col_h:col_h+N-1);
      ## eval fitness.
      [yy,discard_flag ] = problem( N, d, lambda,h );
      if ( discard_flag )
	  %remove this population member.
	  population( idx,col_fit) = -1;
          continue
      end
      fit = fitness( yy );
      population( idx, col_h:col_h+N-1 ) = yy.h;
      population( idx,col_fit) = fit;
      population( idx,col_T) = yy.T;
      population( idx,col_TL) = yy.TL;
      population( idx,col_focus) = yy.z0;
 end
 
 ## remove the 'discard_flag population'
 idx = find(population(:,col_fit) == -1 );
 population(idx,:)=[];
 
 ## calc survival probabilities.
 ## probability of survival
 population(:,col_prob) =   population(:,col_fit)./sum( population(:,col_fit) );

 ## terminate least fit / culling
 [val,idx] = sort(population(:,col_prob),'descend'); 
 population(idx(Npop+1:end),:)=[];
 %population

unfit_population = []; #collects ancestors
gen = 1;
while gen < Ngen
	
	 ## re-calc survival probabilities.
	 tot = sum( population(:,col_fit) );
	 prob_survival = population(:,col_fit)./tot;
	 population(:,col_prob) = prob_survival;
        
	## resort
	[val,idx] = sort(population(:,col_prob),'descend');
        fprintf('Fittest member @ gen %g\n',gen);
        val=population(1,:);
        val(1:2)
	
	## clone top-1-fit of this population.
	##cloned_pop = population(idx(1:clone_rate),:);
	cloned_pop = [];
        NewPop = [];

    % iteration to achieve proper legal mutation populations sizes
    while( rows(NewPop) < round(Npop*(xover_prob + mutate_prob)) )
       
	## roulette wheel selection
	Nlarge=100*rows(population);
	prob_survival =  \
		population(:,col_fit)./sum( population(:,col_fit) );
	fitness_seq = weightedsel( prob_survival, Nlarge);
	
	## mutate 20% population members at one location in pool by their p.d.f.
	## We need mutate_prob*rows(population) chosen from rows(population) to 1, uniquely
	mutate_idx = randint(round(rows(population)*(xover_prob+mutate_prob)),1,rows(population),1,true);
	mutate_idx = unique( fitness_seq( mutate_idx ) );
	mutated_pop = population(mutate_idx,:);
        for idx=0:N-1
	     mutated_pop( : , col_d + idx ) = ga_mutate_binvec( mutated_pop( :, col_d + idx ) );
             ## random search allowed.
	     ## mutated_pop( : , col_h + idx ) = mutated_pop( :, col_h + idx ) + (1-2*(rand()>0.5))*randint(1,1,10,1);
	     ## mutated_pop( : , col_h + idx ) = ga_mutate_binvec( mutated_pop( :, col_h + idx ) );
	     mutated_pop( : , col_h + idx ) =  mutated_pop( :, col_h + idx ) + (1 - 2*(rand()>0.5));
        end
	mutated_pop( : , col_lambda ) = ga_mutate_binvec( mutated_pop( :, col_lambda ) );

	
    ## cross-over same members as mutated
    xover_pop = mutated_pop;
    for idx=0:N-1
	xover_pop( : , col_d + idx ) = ga_xover_binvec( xover_pop( :, col_d + idx ) );
	xover_pop( : , col_h + idx ) =  xover_pop(:,col_h+idx) + \
		bitxor( xover_pop( :, col_h + idx )*(1 - 2*(rand()>0.5)), 7 );
    end
    xover_pop( : , col_lambda ) = ga_xover_binvec( xover_pop( :, col_lambda ) );
    
    ## evaluate fitness
    new_pop = [ xover_pop; cloned_pop];

    ## remove n =0 case, a known minimum
    ## remove illegal data point solutions
    ## remove (lambda, d) < 200 as well as h < 30.
    idx = find( new_pop(:,col_lambda) < 200 );
    new_pop(idx,:) = []; 

    rmlist = zeros(rows(new_pop),1);
    for idx = col_d:col_d+N-1
      idy = find( new_pop(:,idx) < 10 );
      rmlist(idy) = 1;
    end
    new_pop(find(rmlist>0),:) = [];
    
    rmlist = zeros(rows(new_pop),1);
    for idx = col_h:col_h+N-1
      idy = find( new_pop(:,idx) < 10 );
      rmlist(idy) = 1;
    end
    new_pop(find(rmlist>0),:) = [];

    for idx = 1:size(new_pop,1)
          val = new_pop(idx,:);
          d = val(col_d:col_d+N-1); lambda = val(col_lambda);  h = val(col_h:col_h+N-1);
          # eval fitness.
	  [yy,discard_flag ] = problem( N, d, lambda,h ); 
          if ( discard_flag )
	     %remove this population member.
	     new_pop( idx,col_fit) = -1; 
             continue
          end
	  fit = fitness( yy );
          new_pop( idx, col_h:col_h+N-1 ) = yy.h;
          new_pop( idx,col_fit) = fit;
          new_pop( idx,col_T) = yy.T;
          new_pop( idx,col_TL) = yy.TL;
     end
     
    ## remove the 'discard_flag population'
    idx = find(new_pop(:,col_fit) == -1 );
    new_pop(idx,:)=[];
    NewPop = [NewPop; new_pop];
   end
   
    ## consolidate population
    population = [population; NewPop ];
    
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
tot = sum( population(:,col_fit) );
prob_survival = population(:,col_fit)./tot;
population(:,col_prob) = prob_survival;

## resort and retain fittest.
[val,idx] = sort(population(:,col_prob),'descend'); 

## topmost is the best now
population = population(idx,:);

## display
fprintf('generation = %d, best-solution\n',gen);
fprintf('%g \n',population(1,:));
end
