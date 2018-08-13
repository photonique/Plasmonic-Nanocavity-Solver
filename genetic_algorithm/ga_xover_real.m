%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%% (C) 2007-2012 Muthiah Annamalai <muthiah.annamalai@uta.edu>
%%
%% 
%% 12/14/07: deal with xover ga binary vectors.
%% given a vec (row of vectors) we xover two vectors
%% at a random xover site. If the ga_vec is a series of numbers
%% then we use binvec() to convert them to a series of 
%% vectors and back! we expect ga_vec to be an even length
%% for xover operations. If length of ga_vec is odd, the last
%% is not xover with anything, and returned as such.
%% 
%% Input : generic GA xover vector - 1xN or may MxN.
%% Output : val - 1xFloor[N/2] or MxFloor[N/2].
%% 
function val = ga_xover_real(ga_vec)

   if ( nargin < 1 )
      print_usage();
   end

   % mection for dealing with vector bin-strings.
   r = rows(ga_vec);
   c = columns(ga_vec); #just need r/2 xover sites.
   %xover_site = min( fix(rand(r/2,1)*c) + 1, c);
   
   %
   %val = ga_vec;
   %
   %for rx = 1:floor(r/2)
   %  v1 = 2*rx -1 ; v2 = 2*rx;
   %  val(v1,:) = ga_mixup_weights(  ga_vec(v1,:), ga_vec(v2,:) );     
   %  val(v2,:) = ga_mixup_weights(  ga_vec(v1,:), ga_vec(v2,:) );
   %end
   %
   
   r2 = floor(r/2);
   val = ga_vec(:,1);
   tmp = ga_mixup_weights( ga_vec(1:r2,:), ga_vec(r2+1:2*r2,:) );
   val(1:r2) = tmp(:,1);
   tmp = ga_mixup_weights( ga_vec(r2+1:2*r2,:), ga_vec(1:r2,:));
   val(r2+1:2*r2) = tmp(:,1);  
   
   return
end

% xup real vector elements in the two parent into an offspring.
% vA, vB have same number of vector elements.
% return just one column ?
function vO = ga_mixup_weights( vA, vB )
        SZ = size(vA);        
        vA = vA(:);
	vB = vB(:);
	choice = [1 + floor(rand(length(vA),1)*2)];	
	choice = sub2ind([length(vA),3],[1:length(vA)]',choice);
	frac = rand(length(vA),1);
	opts = [vA.*(1+frac), (vA.*frac+(1-frac).*vB), vB.*(1-frac)];
	vO  = opts(choice);
	vO = reshape(vO,SZ(1),SZ(2));
   return;
end
