%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
##
## 01/26/09: xover the permutations, applying one to the other.
##           permutations have to be in 1:length(X) format.
## 
function Val = ga_xover_perm( P, Q )

   if ( nargin < 2 )
      print_usage();
   end

   N = columns( P );
   Val = zeros( size( P ) );

   for idx = 1:rows( P )
     p = P( idx, : );   q = Q( idx, : );

     if ( rand() >= 0.5 )
       val = p ( q );
     else
       val = q ( p );
     end     

     Val( idx, : ) = val;
   end

   return;
end
