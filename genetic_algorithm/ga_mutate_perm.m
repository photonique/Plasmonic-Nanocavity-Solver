%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
##
## 02/03/09: Change order of limits, in calls to randint.
##
## 01/26/09: mutate the permutation, swapping elements.
##           permutations have to be in 1:length(X) format.
## 
function val = ga_mutate_perm( P )

   if ( nargin < 1 )
      print_usage();
   end
   
   N = columns( P );
   val = zeros( size( P ) );

   for idx = 1:rows( P )
     Q = randint( 1, 2, N, 1, true );
     val( idx, : ) = P( idx, : ) ;
     val( idx, Q(1:2) ) = val ( idx, Q(2:-1:1) );
   end
   
   return;
end
