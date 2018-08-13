%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% 03/09/09: Correct for 0-groove case by treating separately.
% 12/19/07: solver for aperiodic case
% 12/15/07: convert solver for use with a GA.
% 

function [transmi_total,e_alpha]=ttotal_aperiodic_solver( n, D, lambda, a, h )
   
   if ( nargin < 3 ), error('aperiodic_solver_scaled_nanocavity'), end;
   if ( nargin < 4 ),   a=40,  end;
   if ( nargin < 5 ), h=10.1:2:(lambda/4+5);, end;

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


   if ( n > 0 )
     rhs=zeros(2*n+1,1);
     rhs(n+1)=2*sqrt(-1);
     eps=ones(2*n+1,1)/tan(k*h);
     eps(n+1)=-sqrt(-1);
     eps1=diag(eps);
     e_alpha=inv(g-eps1)*rhs;
   else
     rhs = 2*sqrt(-1);
     eps1 =-sqrt(-1);
     e_alpha=rhs/(g-eps1);
   end
   
   transmi_total=1-(abs(1-e_alpha(n+1,:))).^2; 
   return
end
%%
%% ttotal_aperiodic_solver( 0, 200, 200, 50, 50.1 )
%% 
