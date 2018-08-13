%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
##
## 02/02/09: If isunique is used, a row-vector or column-vector is assumed.
##           A greedy approach to isunique is followed.
## 
## 12/14/07: 
## pick a random integer using rand() within the 
## range [top, bot] (both inclusive). Also
## constraint is top > bot. row & col specifiy
## the shape of the matrix. Unique flag applies
## to the vectors only.
##
function x = randint( row, col, top, bot, isunique, setval)
    if ( nargin < 2 ), print_usage(), end;
    if ( nargin < 5 ),  isunique = 0; , end;
    if ( nargin < 4 ),  bot = 0; , end;
    if ( nargin < 3 ),  top = 1; , end;
    if ( nargin < 6 ),  setval = [ bot:top ]; , end;

    if ( isunique && ( row*col > numel( setval ) ) )
      warning(' Unique number of elements not generated. Pigeon hole  principle disallows this.' );
    end
    
    x = rand( row, col );
    T = max(top,bot); B = min(top,bot);
    x = (T-B)*x + B;
    x = round(x);

    ## FIXME: randint() on unique vector set.
    ## no equivalent exists for the matrices.
    if ( isvector( x ) && isunique )        
        x = unique( x );
        if ( length(x) < row*col )
	    usedup = x;
	    setval = setdiff( setval, usedup );

	    if ( rows(x) == 1 )
              x = [x, randint( row, col - numel( x ), top, bot, \
			      isunique, setval)]; 
	    else
              x = [x; randint( row - numel( x ), col, top, bot, \
			      isunique, setval)]; 
	    end

        end
        return
    end
    return
end
