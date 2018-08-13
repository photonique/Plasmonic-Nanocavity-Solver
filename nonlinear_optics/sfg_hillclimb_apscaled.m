%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C) 2009 Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
##  
## 02/06/09 : Imply the requirement of the Ttotal @ 775 to be < 1e-2, and
##            no more incentive to reach higher up.
## 
## 02/05/09 : Add random number offsets / random walk type operations 
## on the local search operator to generate new designs. Change dd, dh.
## 
##  Sum Frequency Generation ( SFG ) optimization as a two wavelength
##  joint solver.
##  
##  Wp + Ws =  Wo, where 
##  (signal) Ls = 4100nm, and 
##  (pump)  Lpump = 775nm add in the Chi-2 crystal and generate
##  the output (idler) Lout = 650nm.
## 
##      +----------+-------+ 
##      | Ls       ->  Tmax|
##      +------------------+
##      | Lp       ->  Rmax|
##      +----------+-------+
## 
## LambdaSP = signal, and  pump wavelength combinations.
## 
function population = sfg_hillclimb_apscaled( Ngen, N, xover_prob, mutate_prob, LambdaSP, best_solution )

if ( nargin < 5 ), print_usage(), end;
if ( nargin != 6 )
  best_solution = [];
end

fprintf("HillClimb: N = %g, Ngen = %g\n",N,Ngen);
## keep this even #
Npop = 50; UNFIT_MAX = 150; unfit_population = [];
Hpop = 10;

clone_rate = 1;

## wavelengths of signal, pump and idler.
Lsig = LambdaSP( 1 );
Lpump = LambdaSP( 2 );

## 
## fitness corresponds to maximizing signal tx, pump and idler
## confinement upto a 1e-2 size (see table above).
## 
## This is best fitness function IMHO.
fitness_max =  Inf;

fitness = @( Tsig, Ttot_pump, TLpump )( 10*Tsig*5*(1.1 - (Ttot_pump>0.1)) + (Ttot_pump < 0.1)*10*TLpump );

 ## fixed parameter a=40nm is implicit to solver.
 problem_sig = @(N,D,Lambda,H) ( aperiodic_depthsolver_scaled(N,D,Lambda, 40,H) );
## problem_pump =  @(N,D,Lambda,H) ( aperiodic_depthsolver_scaled_ttotal(N,D,Lambda, 40,H) );

 ## keep data format compatible with the OPA-3wavelength searches.
 
 ## various alleles & their bit allocations
 ## both col_T, col_TL contain the data Ttotal.
 col_Tsig = 1;
 col_Tpump = 2;
 col_TLpump = 3;
 col_lambda = 4;
 col_prob = 5; ## 10bit
 col_fit = 6;
 col_d = 7; #col_d = 7:7+n-1 cols are filled 
 col_h = col_d + N; ## 8bit
 col_total = col_h;

 dh = 10; dd = 10;
 
 population = zeros(Npop*2,col_total);
 population(:,col_fit) = zeros( Npop*2, 1);
 population(:,col_d:col_d+N-1) = randint( Npop*2, N , round(Lsig/2) , round(Lpump));
 population(:,col_h) = randint( Npop*2, 1, round(Lpump), round(Lpump/10));

 population = [best_solution; population];

 ## evaluate fitness.
 for idx = 1:rows(population)
      val = population(idx,:);
      d = val(col_d:col_d+N-1); h = val(col_h);

      ## eval fitness.
      [Ysig,flag_s] = problem_sig( N, d, Lsig,repmat(h,[1,N]) );
      [Ypump,flag_p] = problem_sig( N, d, Lpump,repmat(h,[1,N]) );
      
      ## ifonly all the three flags turn to be true:
      if ( ~( flag_p || flag_s ) )
        Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	population( idx,col_fit ) = fitness( Ysig.T,Ttotal_pump, Ypump.TL );
        fprintf("T %g, T %g, TL %g , Fitness %g\n",Ysig.T,Ypump.T,Ypump.TL, \
                population( idx,col_fit ));
	population( idx, col_Tsig ) = Ysig.T;
	population( idx, col_Tpump ) = Ttotal_pump;
	population( idx, col_TLpump ) = Ypump.TL;
      else
	fprintf('Failed strategy!')
	population(idx,col_fit)=-1;
      end
 end

## remove the 'discard_flag population'
idx = find(population(:,col_fit) == -1 );
population(idx,:)=[];

gen = 1;
while ( gen < Ngen )
   
    ## calculate the probabilities  
    population(:,col_fit) =\
	population(:,col_fit)./sum(population(:,col_fit));
    
    ## resort
    [val,idx] = sort(population(:,col_prob),'descend');
    fprintf('Fittest member @ gen %g\n',gen);
    val=population(idx(1),:);
    val
    
    ## the hillclimbing part.
    for idx = 1:min(Hpop,rows(population))
          val = population(idx,:);
          d = val(col_d:col_d+N-1); h = val(col_h);
	  ## generate 2N variations of the params.
	  ## params = [h,d+dd; h,d-dd; h - dh, d; h + dh, d; zeros(2*N,1+N)];
	  params = [ h, d + randint( 1, 1,dd,-dd);
		    h, d + randint( 1, 1,dd,-dd);
		    h - randint(1,1,dh,-dh), d;
		    h + randint(1,1,dh,-dh), d;
		    zeros(2*N,1+N) ];
	  
	  ## sum & difference of the lengths
	  for idy = 1:N
	    dtemp = d; dtemp(idy) = dtemp(idy) - randint( 1, 1, 1, dd);
	    htemp = h + randint(1,1,dh,-dh);
   	    params(2+2*idy-1,:)=[h,dtemp];
   	    params(2+2*idy,:)=[htemp,d];
	  end
	  
	  for idy = 1:rows(params)
	    ## eval fitness.
	    val = params(idy,:);
 	    h = val(1); d = val(2:end);
	    
	    ## remove edge cases from parameter.
	    if ( min(h) < 10 || min(d) < 50 )
	         continue;
	    end
	    
	    ## eval fitness.
	    [Ysig,flag_s] = problem_sig( N, d, Lsig,repmat(h,[1,N]) );
	    [Ypump,flag_p] = problem_sig( N, d, Lpump,repmat(h,[1,N]) );
	    
	    ## ifonly all the three flags turn to be true:
	    if (  ~( flag_p || flag_s )  )
                Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
                fit = fitness( Ysig.T, Ttotal_pump, Ypump.TL );
	        
                if ( fit > population(idx,col_fit) )
	             
	             ## terminate least fit / culling 
	             unfit_population = [unfit_population; population(idx,:)];
	             ## replace current population with fittest
	             population(idx,col_fit) = fit;
	             population(idx,col_d:col_d+N-1) = d;
	             population(idx,col_h) = h;
	             population( idx,col_Tsig) = Ysig.T;
                     population( idx, col_Tpump ) = Ttotal_pump;
 	             population( idx, col_TLpump ) = Ypump.TL;
	             
 	         end
	      end
            end
	 end
	 
	 ## genetic algorithm section
	 ## the genetic algorithm part
	 ## re-calc survival probabilities. 
	 population(:,col_prob) = population(:,col_fit)./sum( population(:,col_fit) );
	
	## clone top-1-fit of this population.
	cloned_pop = population(idx(1:clone_rate),:);
	
	## roulette wheel selection
        Nlarge=100*rows(population);
	prob_survival = population(:,col_fit)./sum( population(:,col_fit) );
	
	fitness_seq = weightedsel( prob_survival, Nlarge);
	
	## mutate 20% population members at one location in pool.
	## by their p.d.f
	mutate_idx = randint(round(Npop*mutate_prob),1,Nlarge,1,true);
	mutate_idx = unique( fitness_seq( mutate_idx ) );
	
	mutated_pop = population(mutate_idx,:);
	
    for idx=0:N-1
	  mutated_pop( : , col_d + idx ) = ga_mutate_binvec( mutated_pop( :, col_d + idx ) );
    end
        
	mutated_pop( : , col_h ) = ga_mutate_binvec( mutated_pop( :, col_h ) );
	
	xover_idx = randint(round(Npop*xover_prob),1,Nlarge,1,true);
	xover_idx = unique( fitness_seq( xover_idx ) );
	
	xover_pop = population(xover_idx,:);
	
        
    for idx=0:N-1
    	xover_pop( : , col_d + idx ) = ga_xover_binvec( xover_pop( :, col_d + idx ) );
    end
        
	xover_pop( : , col_h ) = ga_xover_binvec( xover_pop( :, col_h ) );
        
	## evaluate fitness of xover+mutated and cloned population
	new_pop = [ xover_pop; mutated_pop; cloned_pop];

	## remove n =0 case, a known minimum
	## remove illegal data point solutions
	## remove (lambda, d) < 200 as well as h < 30.
	idx = find( new_pop(:,col_h) < 10 );
	new_pop(idx,:) = [];

	for idx = col_d:col_d+N-1
            idy = find( new_pop(:,idx) < 50 );
	    new_pop(idy,:) = [];
	end

	## evaluate the clone population
        for idx = 1:size(new_pop,1)

          val = new_pop(idx,:);
          d = val(col_d:col_d+N-1);  h = val(col_h);

	    ## eval fitness.
	    [Ysig,flag_s] = problem_sig( N, d, Lsig,repmat(h,[1,N]) );
	    [Ypump,flag_p] = problem_sig( N, d, Lpump,repmat(h,[1,N]) );
	    
	    ## ifonly all the three flags turn to be true:
	    if (  ~( flag_p  || flag_s )  )
                Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
                fit = fitness( Ysig.T,Ttotal_pump, Ypump.TL );
	        ## replace current population with fittest
	        new_pop(idx,col_fit) = fit;
	        new_pop( idx,col_Tsig) = Ysig.T;
                new_pop( idx, col_Tpump ) = Ttotal_pump;
                new_pop( idx, col_TLpump ) = Ypump.TL;
            else
		fprintf('Failed strategy!')
		new_pop(idx,col_fit)=-1;
	    end
        end

        ## remove the 'discard_flag population'
        idx = find(new_pop(:,col_fit) == -1 );
        new_pop(idx,:)=[];

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
	     [val,idx] = sort(unfit_population(:,col_fit),'descend');
	     unfit_population = unfit_population(idx(1:UNFIT_MAX),:);
	end

    ## remove duplicates in the population and fill with topmost 
    ## unused populations.
    ## size(unfit_population)
    ## size(population)
    [upop,idpop,idupop]=unique(population,'rows');
    population = upop;

    if ( Npop > rows(population) )
	idxmissing = Npop - rows(population)
	if (idxmissing < rows(unfit_population) )
           population = [population; unfit_population(1:idxmissing,:)];
           unfit_population=unfit_population(idxmissing+1:end,:);
	end
    end
    
    [val,idx] = sort(population(:,col_fit),'descend');
    population = population(idx,:);
    if( rows(population) > Npop )
      unfit_population = [unfit_population; population(Npop+1,:)];
      population = population(1:Npop,:);
    end
    
    ## termination condition
    if ( gen > Ngen/2 && max(population(:,col_fit)) > fitness_max )
	break;
    end
    
    ## keep going
    gen = gen + 1;
end

## resort and retain fittest.
[val,idx] = sort(population(:,col_fit),'descend'); 

## topmost is the best now
population = population(idx,:);
population = population(1:min(Npop,rows(population)),:);

## display
fprintf('generation = %d, best-solution\n',gen);
fprintf('%g \n',population(1,:));
end
%%
%% sfg_hillclimb_apscaled(Ngen,N,Lambda,best_solution)
%% sfg_hillclimb_apscaled(2,2,0.2,0.1,[4100 775])
%% 


