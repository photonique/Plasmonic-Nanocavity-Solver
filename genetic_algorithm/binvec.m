%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%
% (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
% University of Texas, Arlington
%
%
% given a decimal vector, converts to the
% binary vector. This is a matrix with each row
% being the 'vector' form of the corresponding 
% decimal value. Expects to work on +ve integers.
%
% sum along the rows, we get back same numbers;
% example: sum((binvec(0:7).*[repmat([4 2 1],[8,1])]),2) 
%

% FIXME: handle negative numbers as well.
function bvec=binvec(dec_vec)
     maxlen=ceil(log2(max(dec_vec)+1));
     x=[]; bvec=zeros(length(dec_vec),maxlen);
     for idx=maxlen:-1:1
         tmp=mod(dec_vec,2);
         bvec(:,idx)=tmp.';
         dec_vec=(dec_vec-tmp)./2;
     end
     return
end
