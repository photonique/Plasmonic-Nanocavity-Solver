%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## (C) 2008, 2009 Muthiah Annamalai <muthiah.annamalai@uta.edu>
## 
## 03/03/09: Raise aspiration levels to 20% each. Change dd, dh.
##           Raise spread of d and h variable. Reduce back a little bit.
##           as no improvement was found.
##           
## 02/27/09: Modify the hillclimbing use of randint.
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
## 12/08/08: New fitness to add another product function, hopefully
##           newer results.
## 
## 12/07/08: New choice of fitness functions based on interpolation points.
## 
## 11/27/08: Adding the forgotten GA searches, a new low.
## 
## 11/15/08: solver for OPA: (signal gets amplified by 
##           a high frequency pump, and idler and pump
##           need to have a large reflection coefficient).
##  
##  Wp = Wo + Ws, where 
##  (idler) Lout = 950nm, (signal) Ls = 4100nm, and 
##  (pump)  Lpump = 775nm
##  
##      +----------+-------+ 
##      | Ls       ->  Tmax| 
##      +------------------+
##      | Lp & Lo  ->  Rmax|
##      +----------+-------+
## 
## LambdaSPO = signal, pump and output/idler wavelength combinations.
## 
function population = \
      opa_hillclimb_lambda_aperiodic_scaled(Ngen,N,xover_prob, \
					    mutate_prob, LambdaSPO,best_solution)
## FIXME: Add similarity removal in selection-culling process.

if ( nargin < 5 ), print_usage(), end;
if ( nargin != 6 )
  best_solution = [];
end

fprintf("HillClimb: N = %g, Ngen = %g\n",N,Ngen);

## keep this even #
Npop =100; UNFIT_MAX = 50; unfit_population = [];
Hpop =5;

clone_rate = 1;

## wavelengths of signal, pump and idler.
Lsig = LambdaSPO(1); Lpump = LambdaSPO(2); Lopidler = LambdaSPO(3);

## 
## fitness corresponds to maximizing signal and idler tx to gaussian mode, 
## pump and confinement for Ttotal -> 1e-2 size, at the highest T/L.
## 
fitness_max =  Inf;
# fitness = @( Tsig, Ttot_pump, \
# 	    TLpump, Tidler )( 50*Tidler*(1.1-(Ttot_pump>0.1)) \
#			     +50*Tsig + (Ttot_pump < 0.1)*10*TLpump );
 fitness = @( Tsig, Ttot_pump, \
 	    TLpump, Tidler )( 1/(1e-3 + norm( [Tidler,Tsig] - [0.15,0.15],2)) ); #\
 			     #+0*(Ttot_pump < 0.1)*10*TLpump );
#
# ((Tidler + 0.1)*Tsig)/(Tidler + Tsig) 
# 50*(Tidler + Tsig)/(1e-3 + Tidler + Tsig) ; fitness 3 really bad
# 50*(Tidler + Tsig) ; fitness 2
# fitness = @( Tsig, Ttot_pump, \
# 	    TLpump, Tidler )( 10*(Tidler + Tsig + 0.25*TLpump)*(1.1 - (Ttot_pump > 0.1)));
#
 
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
 col_d = 7; #col_d = 7:7+n-1 cols are filled 
 col_h = col_d + N; ## 8bit
 col_total = col_h;

 dh = 10; dd = 10;% or 3, 5 or 5, 10 pair.
 
 population = zeros(Npop*2,col_total);
 population(:,col_fit) = zeros( Npop*2, 1);
 population(:,col_d:col_d+N-1) = randint( Npop*2, N , round(Lsig/2) , round(Lpump/2));
 population(:,col_h) = randint( Npop*2, 1, round(Lpump/2), round(Lpump/5));
 
 population = [best_solution; population];

 ## evaluate fitness.
 for idx = 1:rows(population)
      val = population(idx,:);
      d = val(col_d:col_d+N-1); h = val(col_h);

      ## eval fitness.
      [Ysig,flag_s] = problem_sig( N, d, Lsig,repmat(h,[1,N]) );
      [Ypump,flag_p] = problem_sig( N, d, Lpump,repmat(h,[1,N]) );
      [Yopidler,flag_op] = problem_sig( N, d, Lopidler,repmat(h,[1,N]) );
      
      ## ifonly all the three flags turn to be true:
      if ( ~( flag_p || flag_op || flag_s ) )
	Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	population( idx,col_fit ) = fitness( Ysig.T,Ttotal_pump, Ypump.TL,Yopidler.T );
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
          d = val(col_d:col_d+N-1); h = val(col_h);
	  ## generate 2N variations of the params.
	  ## params = [h,d+dd; h,d-dd; h - dh, d; h + dh, d; zeros(2*N,1+N)];
	  params = [ h, d + randint( 1, N,dd,-dd);
		    h, d + randint( 1, N,dd,-dd);
		    h - randint(1,1,dh,-dh), d;
		    h + randint(1,1,dh,-dh), d;
		    zeros(2*N,1+N) ];
	  
	  ## sum & difference of the lengths
	  for idy = 1:N
	    dtemp = d; dtemp(idy) = dtemp(idy) + randint( 1, 1, dd, -dd);
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
	    [Yopidler,flag_op] = problem_sig( N, d, Lopidler,repmat(h,[1,N]) );
	    
	    ## ifonly all the three flags turn to be true:
	    if (  ~( flag_p || flag_op || flag_s )  )
 	        Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	        fit = fitness( Ysig.T,Ttotal_pump, Ypump.TL,Yopidler.T );
	        
                if ( fit > population(idx,col_fit) )
	             
	             ## terminate least fit / culling 
	             unfit_population = [unfit_population; population(idx,:)];
	             ## replace current population with fittest
	             population(idx,col_fit) = fit;
	             population(idx,col_d:col_d+N-1) = d;
	             population(idx,col_h) = h;
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
          [Yopidler,flag_op] = problem_sig( N, d, Lopidler,repmat(h,[1,N]) );
      
	  ## ifonly all the three flags turn to be true:
	  if ( ~( flag_p || flag_op || flag_s ) )
	     Ttotal_pump = Ypump.T.*(1 + 1./Ypump.TL);
	     new_pop( idx,col_fit ) = fitness( Ysig.T,Ttotal_pump, Ypump.TL,Yopidler.T );
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
