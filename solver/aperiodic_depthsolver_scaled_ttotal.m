%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% 11/13/08 : Calcualate the minimum transmission for the case 
% 

function [ttotal,discard_flag]=aperiodic_depthsolver_scaled_ttotal( n, D, lambda, a, hrange_cartprod )
   
   if ( nargin < 3 ), error('aperiodic_solver_scaled_nanocavity'), end;
   if ( nargin < 4 ),   a=40,  end;
   if ( nargin < 5 ), h=10.1:2:(lambda/4+5);, end;


   discard_flag = false;
   a0=1400*lambda/(560);%scaled waist of Gaussian mode 
   spot=a0/a; 
   k=2*pi/lambda; 
   kappa=k*a; 
   ds = [-fliplr(cumsum(D)) 0 cumsum(D)];
    
   %The following lines calculate g=g_alpha_beta in spatial domain (very
   %efficient)
    
   p = [0];
   for x = 2:2*n+1
    for y = 1:x-1
     val = (ds(x) - ds(y))/a;
     p = [p; val];
    end
   end
   
   x_step=2/50; %denominator (i.e. number of steps) must be even 
   x=-1:x_step:1; 
   grid_x=ones(size(x)); 
   grid_x(1)=1/2; 
   grid_x(length(x))=1/2; 
   zz=abs(ones(size(p))*x+p*ones(size(x))); 
   zz(1,(length(x)+1)/2)=1; 
   hh=besselh(0,1,kappa*zz); 
   hh(1,:)=hh(1,:)-sqrt(-1)*2*log(zz(1,:))/pi; 
   hh(1,(length(x)+1)/2)=1+sqrt(-1)*2*(log(kappa/2)+0.57721566490153286061)/pi;

   gp=x_step*hh*(grid_x.*(1-abs(x)))';
   gp(1)=gp(1)-3*sqrt(-1)/pi; 
   gp=gp*sqrt(-1)*kappa/2; 
   
   % for certain "bad" parameters model doesnt work at all
   if ( isnan(max(gp)) )
        ttotal = [];
	discard_flag=true;
	 n
	 D
	 lambda
	 hrange_cartprod
	fprintf('Failed!\n')
	return
   end
   
   itr = 1;
   for x = 1:2*n+1
     for y = 1:x
      if ( x == y )
        g(x,y) = gp(1);
      else
        val = gp(itr+1);
        itr = itr+1;
        g(x,y) = val;
        g(y,x) = val;
      end
     end
   end 

   rhs=zeros(2*n+1,1); 
   rhs(n+1)=2*sqrt(-1); 

   for i=1:rows(hrange_cartprod)
       eps=ones(2*n+1,1)./[tan(k.*fliplr(hrange_cartprod(i,1:n))), 1, tan(k.*hrange_cartprod(i,1:n))].';
       eps(n+1)=-sqrt(-1);
       eps1=diag(eps);
       G = (g-eps1);
       
       if ( max(max(isnan(G)+isinf(G))) >= 1 || cond(G) > 1e3 )
	 ttotal=[];
	 discard_flag=true;
	 n
	 D
	 lambda
	 hrange_cartprod
	 fprintf('Failed!\n')
	 return
       end       
       e_ans =  G\rhs;
       e_alpha(:,i)=e_ans;
   end
   ttotal=1-(abs(1-e_alpha(n+1,:))).^2;
   discard_flag = false;
   return
end
%
% aperiodic_depthsolver_scaled_ttotal(2, [300 300], 324, [40 50] )
%
