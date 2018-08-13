%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 01/29/09: Use the actual function file for rastrigin.
##
## 01/27/09: Add separate populations of Xover & Mutation.
##
##
##10/21/08: Compare my GA code to Mathworks GA Toolbox.
##          by solving the canonical Rastrigins function.
##          I needed to represent a "real" parameter space
##          using a discrete representations. Works well!
##
## Following lines visualize the search space:
## [xx,yy]=meshgrid(-2:0.1:2,-2:0.1:2);
## mesh(20+xx.^2 + yy.^2 - 10.*(cos(2*pi*xx) + cos(2*pi*yy)))
##

clear all
close all
more off

## find the global minima (0,0) of the rastrigins function:
## ras = @( x ) (20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) +    cos(2*pi*x(2)) ) );
fitness = @(result) 45-result; ##seek minima of rastrigins function
fitness = @(result) 1/(1+result); ##seek minima of rastrigins function

Ngen=2*25; xover_prob = 0.15;
Npop =50; mutate_prob = 0.1; 
clone_rate = 1; UNFIT_MAX = 1e3;

#record fitness, best solutions as function of generations
Fitness = []; 
BestSol = [];

# use search space as [-2,-2] to [+2,+2] rectangular area

## fixed parameter a=40nm is implicit to solver.
problem = @(x1,x2) ( rastrigin([x1,x2]) );

 ## various alleles & their bit allocations
 col_x1 = 1;
 col_x2 = 2;
 col_fit = 3;
 col_prob = 4;
 col_total = col_prob;

 population = zeros(Npop,col_total);
 ## we have a 'real' space in -2,2, but use 400 no's to represent this.
 scale_x = 50; scale_y = scale_x;
 population(:,col_x1) = randint( Npop, 1, 2*scale_x,-2*scale_x);
 population(:,col_x2) = randint( Npop, 1, 2*scale_x,-2*scale_x);
 ## evaluate fitness.
 for idx = 1:rows(population)
      ##analp(population)
      val = population(idx,:);
      x1 = val(1); x2=val(2);
      # eval fitness.
      yy = problem( x1/scale_x,x2/scale_y); fit = fitness( yy );
      population( idx,col_fit) = fit;
 end
 
 ## calc survival probabilities.
 ## probability of survival
 population(:,col_prob) =   population(:,col_fit)./sum( population(:,col_fit) );

 ## terminate least fit / culling 
 [val,idx] = sort(population(:,col_prob),'descend');
 population(idx(Npop+1:end),:)=[];

 unfit_population = []; #collects ancestors
 gen = 1;
 while gen <= Ngen
	
	 ## re-calc survival probabilities.  
	 population(:,col_prob) = population(:,col_fit)./sum( population(:,col_fit) );

	## resort
	[val,idx] = sort(population(:,col_prob),'descend'); 
    fprintf('Fittest member @ gen %g\n',gen);
	val=population(1,:);
	[val(1)/scale_x, val(2)/scale_y , val(3:end)]
    
	## clone top-1-fit of this population.
	cloned_pop = population(idx(1:clone_rate),:);

	## roulette wheel selection
	Nlarge=100*rows(population);
	prob_survival =  \
		population(:,col_fit)./sum( population(:,col_fit) );

	fitness_seq = weightedsel( prob_survival, Nlarge);
	
	## mutate population members at one location in pool.
	## by their p.d.f
	mutate_idx = randint(round(Npop*mutate_prob),1,Nlarge,1,true);
	mutate_idx = unique( fitness_seq( mutate_idx ) );

	mutated_pop = population(mutate_idx,:);
	mutated_pop( : , col_x1 ) = ga_mutate_binvec( mutated_pop( :, col_x1 ) );
	mutated_pop( : , col_x2 ) = ga_mutate_binvec( mutated_pop( :, col_x2 ) );

	xover_idx = randint(round(Npop*xover_prob),1,Nlarge,1,true);
	xover_idx = unique( fitness_seq( xover_idx ) );

        xover_pop = population( xover_idx,:);
	xover_pop( : , col_x1 ) = ga_xover_binvec( xover_pop( :, col_x1 ) );
	xover_pop( : , col_x2 ) = ga_xover_binvec( xover_pop( :, col_x2 ) );
    
	## evaluate fitness of xover+mutated and cloned population
	new_pop = [ xover_pop; mutated_pop; cloned_pop];


    for idx = 1:size(new_pop,1)
          val = new_pop(idx,:);
          x1 = val(1); x2=val(2);
          # eval fitness.
          yy = problem( x1/scale_x,x2/scale_y); fit = fitness( yy );          
          new_pop( idx,col_fit) = fit;
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
    Fitness = [Fitness;population(1,col_fit)];
    BestSol = [BestSol; population(1,col_x1)/scale_x, population(1,col_x2)/scale_y];
    gen = gen + 1;
 end

 ## re-calc survival probabilities.
 prob_survival = population(:,col_fit)./sum( population(:,col_fit) );
 population(:,col_prob) = prob_survival;

 ## resort and retain fittest.
 [val,idx] = sort(population(:,col_prob),'descend');

 ## scale the parameters.
 population(:,col_x1) = population(:,col_x1)/scale_x;
 population(:,col_x2) = population(:,col_x2)/scale_y;

 ## topmost is the best now
 population = population(idx,:);
 
 ## display
 fprintf('generation = %d, best-solution\n',gen);
 fprintf('%g \n',population(1,:));

subplot(2,1,1)
plot(1:Ngen,Fitness,'-o;fitness vs n#generations;')
subplot(2,1,2)
plot(1:Ngen,BestSol(:,1),'-o;X1 vs n#generations;',1:Ngen,BestSol(:,2),'-o;X2 vs n#generations;')

figure
[xx,yy]=meshgrid(-2:0.1:2,-2:0.1:2);
surf(20+xx.^2 + yy.^2 - 10.*(cos(2*pi*xx) + cos(2*pi*yy)))
