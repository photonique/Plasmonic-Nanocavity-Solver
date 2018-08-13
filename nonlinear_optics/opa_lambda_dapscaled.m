%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C) 2008, 2009 Muthiah Annamalai <muthiah.annamalai@uta.edu>
##
## 04/24/09: Check if all 3 wavelengths pass the validity checks.
## 
## 04/09/09: Use a 2-norm fitness function.
## 
## 04/08/09: Run a search with metrics as
##           Ecav/Ein for both signal 4100 & pump 775. Idler at 955 will
##           just have T (gaussian coupling) as metric. Change fitness fn.
## 
## 04/02/09: Compromise the peak finding algorithm between 3 objectives
##           T(4100,955), and Ecav/Ein.
## 
## 04/01/09: Ecav/Ein = t/(1 - r) = sqrt(Tgauss)/(1 - sqrt(1 - Ttotal));
##
## 03/31/09: Actual correct F-P cavity enhancement of Ecav/Ein, after its
##           verified by Prof. V.
## 
## 03/30/09: Fitness function based new derivation of FP-cavity effect.
##
## 03/27/09: Add a fitness function using Fabry-Perot cavity effect
##           for the values of the 775nm pump max Ecav/Ein inside cavity.
## 
## 03/18/09: Issue being T955 is neglected, while T/L775 works well.
##           So modify fitness function be reducing the scaling.
## 
## 03/17/09: Add all 4 constraints to the OPA problem.
## 
## 03/10/09: Add the T/L from pump wavelength.
## 
## 03/08/09: Reduce the aspiration levels. 
##
## 03/06/09: GA reached aspiration levels to raising it upto 45%.
##
## 03/04/09: Remove repmat on h; h is already 1xN sized.
##
## 03/03/09: Raise aspiration levels to 25% each. reduce dd, dh and Npop.
##           Raise range of d, and h depth variable for a larger spread.
##           Pullback a bit.
## 
## 02/27/09: Complete depth variable search with local search options.
## 02/25/09: Same as older OPA search, but account for depth variations.
## 
## 02/12/09: Search for 3-wavelength, 4constriants; seek maximum T
##           for the signal & idler wavelength. Problem is modified into 
##           the following :
##           
##           Signal = 4100nm. Pump = 775nm. Idler = 650nm.
##           
##           Signal + Idler => Tmax
##           Pump => Rmax ( Ttotal -> 0 ) & T/L max
##           
##           Remove traces of the polyfit() extracted fitness functions.
##           Add random operators in the hillclimbing local search.
##           
##           
function population = \
      opa_lambda_dapscaled(Ngen,N,xover_prob, \
					    mutate_prob, LambdaSPO,best_solution)
## FIXME: Add similarity removal in selection-culling process.

if ( nargin < 5 ), print_usage(), end;
if ( nargin != 6 )
  best_solution = [];
end

fprintf("HillClimb: N = %g, Ngen = %g\n",N,Ngen);

## keep this even #
Npop =50; UNFIT_MAX = 50; unfit_population = [];
Hpop =10; dh = 3; dd = 10;% or 3, 5 or 5, 10 pair.

clone_rate = 1;
aspiration_levels = [ 0.18,0.14];
## wavelengths of signal, pump and idler.
Lsig = LambdaSPO(1); Lpump = LambdaSPO(2); Lopidler = LambdaSPO(3);

## 
## fitness corresponds to maximizing signal and idler tx to gaussian mode,
## pump and confinement for Ttotal -> 1e-2 size, at the highest T/L.
## 
fitness_max =  Inf;

##fitness = @( Tsig, Tpump, Tidler )  \
##    ( 1/(1e-3 + norm( [Tsig.Ecav_Ein,Tpump.Ecav_Ein,Tidler.T] - [50,50,0.15],2) ) );

## fitness = @( Tsig, Tpump, Tidler )  \
##    ( (Tsig.Ecav_Ein > 50 ) + (Tpump.Ecav_Ein > 50 ) + ( Tidler.T > 0.1 ));

## tnorm = @( val, targ ) ( 1/norm( [val - targ],2 ) );

## fitness =  @( Tsig, Tpump, Tidler )  \
##    ( tnorm(Tsig.Ecav_Ein, 25 ) + tnorm(Tpump.Ecav_Ein, 25 ) + tnorm( Tidler.T,0.2));

## fitness = @( Tsig, Tpump, Tidler ) \
##  ( Tsig.T + Tidler.T + Tpump.T + Tpump.Ecav_Ein/10 + Tsig.Ecav_Ein/10 );

fitness = @( Tsig, Tpump, Tidler ) \
  ( Tidler.T + Tpump.Ecav_Ein/10 + Tsig.Ecav_Ein/10 );

## fixed parameter a=40nm is implicit to solver.
 problem_sig = @(N,D,Lambda,H) ( aperiodic_depthsolver_scaled(N,D,Lambda, 40,H) );

 ## various alleles & their bit allocations
 ## both col_T, col_TL contain the data Ttotal.
 col_Tsig = 1;
 col_Tpump = 2;
 col_Topidler = 3;
 col_lambda = 4;
 col_prob = 5; ## 10bit
 col_fit = 6;
 col_d = 7; ## col_d = 7:7+n-1 cols are filled 
 col_h = col_d + N; ## 8bit
 col_total = col_h + N -1; 
 
 population = zeros(Npop*2,col_total);
 population(:,col_fit) = zeros( Npop*2, 1);
 population(:,col_d:col_d+N-1) = randint( Npop*2, N , round(Lsig/2) , round(Lpump/2));
 population(:,col_h:col_h+N-1) = randint( Npop*2, N, round(Lsig/2), round(Lpump/10));
 
 population = [best_solution; population];
 
 ## evaluate fitness.
 for idx = 1:rows(population)
      val = population(idx,:);
      d = val(col_d:col_d+N-1); h = val(col_h:col_h+N-1);
      
      ## eval fitness.
      [Ysig,flag_s] = problem_sig( N, d, Lsig,h );
      [Ypump,flag_p] = problem_sig( N, d, Lpump,h );
      [Yopidler,flag_op] = problem_sig( N, d, Lopidler,h );
      
      ## ifonly all the three flags turn to be true:
      if ( ~( flag_p || flag_op || flag_s ) )
	Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	population( idx,col_fit ) = fitness( Ysig,Ypump,Yopidler );
	population( idx,col_Tsig) = Ysig.T;
	population( idx,col_Tpump) = Ttotal_pump;
	population( idx,col_Topidler) = Yopidler.T;
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
    population(:,col_fit) = \
	population(:,col_fit)./sum(population(:,col_fit));
    
    ## resort
    [val,idx] = sort(population(:,col_prob),'descend');
    fprintf('Fittest member @ gen %g\n',gen);
    val=population(idx(1),:);
    val

    ## the hillclimbing part.
    for idx = 1:min(Hpop,rows(population))
          val = population(idx,:);
          d = val(col_d:col_d+N-1); h = val(col_h:col_h+N-1);
	  ## generate 2N variations of the params.
	  ## params = [h,d+dd; h,d-dd; h - dh, d; h + dh, d; zeros(2*N,1+N)];
	  params = [ h, d + randint( 1, N,dd,-dd);
		    h, d + randint( 1, N,dd,-dd);
		    h - randint(1,N,dh,-dh), d;
		    h + randint(1,N,dh,-dh), d;
		    zeros(2*N,N+N) ];
	  
	  ## sum & difference of the lengths
	  for idy = 1:N
	    dtemp = d; dtemp(idy) = dtemp(idy) + randint( 1, 1, dd, -dd);
	    htemp = h + randint(1,N,dh,-dh);
   	    params(2+2*idy-1,:)=[h,dtemp];
   	    params(2+2*idy,:)=[htemp,d];
	  end
	  
	  for idy = 1:rows(params)
	    ## eval fitness.
	    val = params(idy,:);
 	    h = val(1:N); d = val(N+1:end);
	    
	    ## remove edge cases from parameter.
	    if ( min(h) < 10 || min(d) < 50 )
	         continue;
	    end


	    ## eval fitness.
	    [Ysig,flag_s] = problem_sig( N, d, Lsig,h );
	    [Ypump,flag_p] = problem_sig( N, d, Lpump,h );
	    [Yopidler,flag_op] = problem_sig( N, d, Lopidler,h );
	    
	    ## ifonly all the three flags turn to be true:
	    if (  ~( flag_p || flag_op || flag_s )  )
 	        Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	        fit = fitness( Ysig,Ypump,Yopidler );
	        
                if ( fit > population(idx,col_fit) )
	             
	             ## terminate least fit / culling 
	             unfit_population = [unfit_population; population(idx,:)];
	             ## replace current population with fittest
	             population(idx,col_fit) = fit;
	             population(idx,col_d:col_d+N-1) = d;
	             population(idx,col_h:col_h+N-1) = h;
	             population( idx,col_Tsig) = Ysig.T;
	             population( idx,col_Tpump) = Ttotal_pump;
	             population( idx,col_Topidler) = Yopidler.T;
	             
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
     mutated_pop( : , col_d + idx ) = ga_mutate_binvec( mutated_pop( :, col_d + idx));
     mutated_pop( : , col_h + idx ) = ga_mutate_binvec( mutated_pop( :, col_h + idx));
    end
        
	
	xover_idx = randint(round(Npop*xover_prob),1,Nlarge,1,true);
	xover_idx = unique( fitness_seq( xover_idx ) );
	
	xover_pop = population(xover_idx,:);
	
        
    for idx=0:N-1
    	xover_pop( : , col_d + idx ) = ga_xover_binvec( xover_pop( :, col_d + idx ) );
	xover_pop( : , col_h + idx ) = ga_xover_binvec( xover_pop( :, col_h + idx ) );
    end
        
	## evaluate fitness of xover+mutated and cloned population
	new_pop = [ xover_pop; mutated_pop; cloned_pop];

	## remove n =0 case, a known minimum
	## remove illegal data point solutions
	## remove (lambda, d) < 200 as well as h < 30.	
	
	for idx = col_d:col_d+N-1
            idy = find( new_pop(:,idx) < 50 );
	    new_pop(idy,:) = [];
	end

	for idx = col_h:col_h+N-1
            idy = find( new_pop(:,idx) < 10 );
	    new_pop(idy,:) = [];
	end
	
	## evaluate the clone population
        for idx = 1:size(new_pop,1)

          val = new_pop(idx,:);
          d = val(col_d:col_d+N-1);  h = val(col_h:col_h+N-1);

          ## eval fitness.
          [Ysig,flag_s] = problem_sig( N, d, Lsig,h );
          [Ypump,flag_p] = problem_sig( N, d, Lpump,h );
          [Yopidler,flag_op] = problem_sig( N, d, Lopidler, h );
      
	  ## ifonly all the three flags turn to be true:
	  if ( ~( flag_p || flag_op || flag_s ) )
	     Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	     new_pop( idx,col_fit ) = fitness( Ysig, Ypump,Yopidler );
  	     new_pop( idx,col_Tsig) = Ysig.T;
	     new_pop( idx,col_Tpump) = Ttotal_pump;
	     new_pop( idx,col_Topidler) = Yopidler.T;
          else
	     fprintf('Failed strategy!')
	     new_pop(idx,col_fit) = -1;
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
%% opa_hillclimb_lambda_aperiodic_scaled(Ngen,N,0.2,0.1,Lambda,best_solution)
%% opa_hillclimb_lambda_aperiodic_scaled(2,2,0.2,0.1,[750 300 650])
%% 
