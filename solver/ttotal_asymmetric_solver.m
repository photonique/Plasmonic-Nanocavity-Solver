%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% 08/11/09 : add depths to be asymmetric as well. Remove default depths
%            and slit-width parameters.
% 
% 08/06/09 : convert to asymmetric groove model, when you have
% 
function [transmi_total,e_alpha]=ttotal_asymmetric_solver( n, D, lambda, a, h )
   
   if ( nargin < 5 ), error('asymmetric_solver_scaled(n,D,lambda,a,h'), end;
   
   if ( length(D) != 2*n || length(h) != 2*n )
     error("Too few period lengths or depths for asymmetric periodic solver");
   end
   
   %a0=1400;
   a0=1400*lambda/(560);%scaled waist of Gaussian mode
   spot=a0/a;
   k=2*pi/lambda;
   kappa=k*a;
   ds = [-fliplr(cumsum(fliplr(D(1:n)))) 0 cumsum(D(n+1:2*n))];
   
   %The following lines calculate g=g_alpha_beta in spatial domain (very efficient)
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
   %The following line is needed to obtain e_alpha for a range of depth h
   %h=100; %to compare with k-domain calculations (single h value), uncomment this line
   for i=1:rows(h),
       eps=ones(2*n+1,1)./[tan(k*h(i,1:n)), 0, tan(k*h(i,n+1:2*n))].';
       eps(n+1)=-sqrt(-1);
       eps1=diag(eps);
       e_alpha(:,i)=(g-eps1)\rhs;
   end
   transmi_total=1-(abs(1-e_alpha(n+1,:))).^2;
end
%% 
%% ttotal_asymmetric_solver( 0, 200, 200, 50, 50.1 )
%% 
