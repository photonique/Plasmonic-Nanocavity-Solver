%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% 03/09/09 : Solver that assumes depth [h], period [d], slit width,
% groove widths [a] are variable, with the only symmetry being preserved.
% 
% Input: 
% 
%     n : number of grooves on either side of the slit
%     D : period of grooves
%     lambda : incident wavelength
%     A : slit & groove widths
%     H : depth of grooves
% 
% expect all of the arguments to be 1 x N column vectors
% 

function [transmi_total,e_alpha]=ttotal_allvar_solver( n, D, lambda, A, H )
   
   if ( nargin < 3 ), error('aperiodic_solver_scaled_nanocavity'), end;
   if ( nargin < 4 ), A=repmat(40,[n+1,1]);, end;
   if ( nargin < 5 ), H=repmat(D/2,[2*n+1,1]);, end;
   
   ds = [-fliplr(cumsum(D)) 0 cumsum(D)];
   H = [fliplr(H), 0, H];
   
   if ( n > 0 )
     A = [fliplr(A(2:end)), A];
   end
   
   k=2*pi/lambda;
   kappa=k*A(n+1);

   %The following lines calculate g=g_alpha_beta in spatial domain (very
   %efficient)
   
   p = [0];
   for x = 2:2*n+1
    for y = 1:x-1
     val = (ds(x) - ds(y))/A(n+1);
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

   groove_factor = A(n+1)./(sqrt(A'*A));
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
   
   g = g.*groove_factor;
   
   rhs=zeros(2*n+1,1);
   rhs(n+1)=2*sqrt(-1); 
   
   
   if ( n > 0 )
     rhs=zeros(2*n+1,1);
     rhs(n+1)=2*sqrt(-1);
     eps=ones(2*n+1,1)./tan(k*H');
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
%% ttotal_allvar_solver( 0, 200, 200, 50, 50.1 )
%% = 0.88810
%% 
%% ttotal_allvar_solver( 1, 200, 200, [50 65], [110] )
%% = 0.89572
