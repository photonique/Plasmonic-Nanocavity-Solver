%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## 
## (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
## 
## Rastrigins function, for function optimization.
##
## find the global minima (0,0) of the rastrigins function:
## ras =  (20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) +    cos(2*pi*x(2)) ) );
function val = rastrigin( X )
  if ( numel( X ) < 2 )
    print_usage( );
  end

  val = 20 + sumsq( X ) - 10*sum( cos( 2*pi*X ) ) + 10*numel(X);
  return
end
