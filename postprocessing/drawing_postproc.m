%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
P={}; Q={};
dirn = "/home/muthu/research/genetic-alg/data-092208/"; 
for N = 1:16
          qname = strcat("GA_APERIODIC_N",num2str(N),"-22-Sep-2008");
          vname = strcat("GA_APERIODIC_N",num2str(N),"_22_Sep_2008");
          load(strcat(dirn,qname,".txt"));
     	  eval(sprintf("p=%s(1:1000,:); q=%s(1001:2000,:);",vname,vname))
          P{N}=p(1,:);
          Q{N}=q(1,:);
end

p = zeros(16,8+16);
q = zeros(16,8+16);
for N = 1:16
  p(N,1:N+8-1) = P{N}(1,:);
  q(N,1:N+8-1) = Q{N}(1,:);
end

close all
col_h = 4;
col_d = 8;
a=40; ofst = [0, 0];
for  N = 1:16
     figure
     d=p(N,col_d:col_d+N-1); h=p(N,col_h);
     ofst=drawcavity(N,d,h,40);
     axis( [ -(N+1)*max(d) +(N+1)*max(d) 0 ofst(2) ] )
     print('-dpng',['geom_','N=',num2str(N),'.png'])
     close
end
break

