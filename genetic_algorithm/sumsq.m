%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% (C) 2012 Muthiah Annamalai
% - What Octave does, MATLAB fails to
function v = sumsq( z )
  v = sum(sum(abs(z(:)).^2))
end
