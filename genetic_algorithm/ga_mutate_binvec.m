%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## ## (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
## 12/14/07: deal with mutating ga binary vectors.
## given a vec (row of vectors) we mutate each vector
## at any one point. If the ga_vec is a series of numbers
## then we use binvec() to convert them to a series of 
## vectors and back!
##
function val = ga_mutate_binvec(ga_vec)

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
   c = columns(ga_vec);
   mutate_pts = min( fix(rand(r,1)*c) + 1, c);

   ##FIXME: parallelize the stuff.
   ##mutate_pts = mutate_pts + fix(rand()*c
   ##ga_vec(mutate_pts) = 1 - ga_vec( mutate_pts );  ## columns

   for rx = 1:r
     ga_vec(rx,mutate_pts(rx)) = 1 - ga_vec(rx,mutate_pts(rx));
   end

   val = binval( ga_vec );
   return;
end

