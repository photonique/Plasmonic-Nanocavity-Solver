%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
## 
## 01/29/09 highly nonlinear function
## 
## 0.5 + ( sin^2 ( sqrt( x^2 + y^2 ) ) - 0.5 )
##       -----------------------------------
##          ( 1.0 + 0.001*(x^2 + y^2))^2
##
function v = schaffer_f6( X )
  nr =  (sin(sqrt(sumsq( X ))).^2 - 0.5);
  dr = (1.0 + 0.001*sumsq( X )).^2;
  v = 0.5 + nr./dr;
  return
end
%!assert(schaffer_f6( [-0.66895545, -3.066477] ) , 0.0097159, 1e-5)
%!assert(schaffer_f6( [-1.4643227, -2.7757885] ) , 0.0097159, 1e-5)
