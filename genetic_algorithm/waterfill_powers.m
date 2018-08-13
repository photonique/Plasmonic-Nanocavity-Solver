%%
%% (C) 2007-2012, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% (C) 2012 Muthiah Annamalai <muthuspost@gmail.com>
% University of Texas - Arlington
% 
% 06/21/12 - equalize complex amplitude coefficients to
% the power, correcting for overshoot/undershoot the
% range P0. i.e sum(abs(B).^2)=P0 irrespective of A overflows or underflows.
% Make complex by choosing a phase exp(1i*theta).
% 
% Input: A - unequalized (w.r.t power equivalence) modes 
%        P0 - pump power 
%        tol - target tolerance for pump powers
% 
% Output: B - equalized mode amplitude upto tolerance w.r.t P0.
% 
function B = waterfill_powers( A, P0, tol )
	if ( nargin < 3), tol = 1e-3; end;
	assert( P0 > 0 && tol > 0 );
	SZ = size(A);
	
	B = (abs(A(:)).^2); %r = rand(size(B));
	B(find( B > P0 )) = 0; %clear the outliers.
	
	pexcess = @(M) (P0 - sum(abs(M)));
	pex = pexcess(B);
	% water-filling algorithm, with successive approximation
	% strategy to arrive at converge
	prev = []; itr =1;
	while ( abs(pex) > tol && itr < numel(A) )
	  [val,idx] = max(B);
	  itr = itr +1;	  	  	  
	  B(idx) = max(abs(pex)/2,val/2);	  
	  pex = pexcess(B);
	end
	% force convergence
	[val,idx]=min(B);
	B(idx) = 0;
	pex = pexcess(B);
	B(idx) = pex;
	B = sqrt(B).*exp(1i*angle(A(:)));
	
	% add a shuffling step - so we iron out any holes
	% or bias in the arrangement steps \\//<)(><)(>\\//
	N = numel(A(:));
	r = rand(1,N*5);
	r2 = mod(floor(r*N),N)+1;	
	B = reshape(B,1,numel(B));
	for itr = r2;
		B = [B(itr+1:end),B(1:itr)];
	end
	B = reshape(B,SZ);	
	%B; sum(abs(B.^2));
	return
end

