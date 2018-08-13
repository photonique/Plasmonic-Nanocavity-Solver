%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
## ## (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
## 12/14/07: Muthiah Annamalai
##
## converts the binary vector bvec to a 
## series of values val. This function is a 
## complement operation to binvec() function.
## Operates only on +ve numbers, and uses concept
## of the binary vector from MSB -> LSB.
##
function val = binval( bvec )

	if ( nargin < 1 )
	   print_usage();
	end

	[r,c] = size(bvec);
        val = sum(bvec.*repmat(2.^[c-1:-1:0],[r,1]),2);
end

