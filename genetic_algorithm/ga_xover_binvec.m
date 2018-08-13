%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## ## (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
## 12/14/07: deal with xover ga binary vectors.
## given a vec (row of vectors) we xover two vectors
## at a random xover site. If the ga_vec is a series of numbers
## then we use binvec() to convert them to a series of 
## vectors and back! we expect ga_vec to be an even length
## for xover operations. If length of ga_vec is odd, the last
## is not xover with anything, and returned as such.
## 
function val = ga_xover_binvec(ga_vec)

   if ( nargin < 1 )
      print_usage();
   end

   ## treat case when ga_vec is NOT a binary vector,
   ## and a vector of numbers.
   if ( isvector(ga_vec) && max( ga_vec ) > 1 )
      ga_vec = binvec( ga_vec );
   end

   ## section for dealing with vector bin-strings.
   r = rows(ga_vec);
   c = columns(ga_vec); #just need r/2 xover sites.
   xover_site = min( fix(rand(r/2,1)*c) + 1, c);

   ##FIXME: parallelize the stuff.

   for rx = 1:r/2 
     v1 = 2*rx -1 ; v2 = 2*rx;
     tmp = ga_vec( v1, : );
     ga_vec( v1, xover_site(rx):end) = ga_vec( v2 ,xover_site(rx):end);
     ga_vec( v2, xover_site(rx):end) = tmp(xover_site(rx):end);
   end
   
   ## FIXME: treat the odd last element, thats not xover
   ## to use the first element as a pair.

   val = binval( ga_vec );
   return;
end

