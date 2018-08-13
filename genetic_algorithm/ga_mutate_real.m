%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%% (C) 2007-2012 Muthiah Annamalai <muthiah.annamalai@uta.edu>
%% 
%% 06/18/12: Reals in mutations
%% 
%% We define real-valued mutation as changing upto +/- 10% of real vector value. 
%% Input: 'ga_vec' - standard real number vector
%% Output: mutated version of input - modified by +/-10% of 'ga_vec'
function val = ga_mutate_real(ga_vec)

   if ( nargin < 1 )
      print_usage();
   end
   
   r = rows(ga_vec);
   c = columns(ga_vec);
   
   ga_vec = ga_vec + (-1).^(rand(r,c)>0.5).*(0.1.*rand(r,c).*ga_vec);
   
   val =  ga_vec;
   return;
end

