%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% (C) 2012 Muthiah Annamalai
% The University of Texas - Arlington.
% 
% 06/21/12: 
% Force power equality in the mode amplitudes.
% 
% 06/15/12:
% Function - ga_psa_crossover
% Combines two parents with structure elements,
% 'weights' - a 2D matrix,
% 'mode_size' - a 2x1 vector or scalar,
% to produce offspring.
% 'P0' - sum square of the mode amplitudes should add up to the power
% 
% Returns as many offspring as largest number of parent passed
% in the arguments. The smaller vector parent has its last element
% duplicated to makeup for the larger parent size.
% 
% Input - parentM, parentD - parent structures with fields, 'weights'
%         and 'mode_size'.
%         'P0' - pump power which is criteria for mode weight equivalence.
% Output - 'offspring' - combined  crossover version of two parent structs
%         such that 'weights' are power equivalent to P0.
% 
function offspring = ga_psa_crossover( parentM, parentD, P0 )
	offspring = struct('weights',{},'mode_size',[]);
        L = [length(parentM),length(parentD)];	

	for itr = 1:max(L)
		itrM = max(itr,L(1));
		itrD = max(itr,L(2));
		offspring(itr).weights = ga_psa_mixup_weights( parentM(itrM).weights, parentD(itrD).weights );
		offspring(itr).mode_size = ga_psa_mixup_weights( parentM(itrM).mode_size, parentD(itrD).mode_size);
		offspring(itr).weights = waterfill_powers(offspring(itr).weights,P0);
	end
	return
end

% mixup real vector elements in the two parent into an offspring.
% vA, vB have same number of vector elements.
function vO = ga_psa_mixup_weights( vA, vB )
	SZ = size(vA);
	%vA = reshape(vA,prod(SZ),1);
	%vB = reshape(vB,prod(SZ),1);
	vA = vA(:);
	vB = vB(:);
	choice = [1 + floor(rand(length(vA),1)*2)];	
	choice = sub2ind([length(vA),3],[1:length(vA)]',choice);
	frac = rand(length(vA),1);
	opts = [vA, (vA.*frac+(1-frac).*vB), vB];
	vO  = opts(choice);
	vO = reshape(vO,SZ(1),SZ(2));
	return
end
