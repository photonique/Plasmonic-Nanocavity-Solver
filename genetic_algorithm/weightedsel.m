%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%%
%% 02/02/09 : Alternative weighted selection mechanism, to simulate the
%%            Roulette wheel, similar to code found in D. E. Goldberg's 
%% landmark textbook "Genetic Algorithms".
%% 
%% Given pdf, and number of points to pick,
%% the function returns the index on pdf, as
%% the selected points. choice is a Nx1 vector,
%% with elements containing integers in range
%% 1 <= choice(x) <= length(pdf).
%% 
%% Input: pdf - Prob.density.function to use with same alphabet size 
%%        N - number of elements to choose from distributed as 'pdf'.
%% Output: choice = 1xN  vector with elements or hist(choice) follows 'pdf'.
%% 
function choice = weightedsel( pdf, N )
  
  if ( nargin < 2 )
    print_usage( )
  end

  choice = zeros( N, 1 );

  cdf = cumsum( pdf );

  wheel_pts = rand( N, 1 );

  for idx  = 1: N
    choice( idx ) = min( find( cdf >= wheel_pts( idx ) ) );
  end

  return
end
