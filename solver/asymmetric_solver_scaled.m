%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
% 
% 03/28/10: Matlab compatible
% 09/23/08: solver for asymmetric-aperiodic case
% 09/23/08: asymmetric case where d[i] != d[-i], and d[i] != d[i+1] with
%           completely aperiodic structure of the metal nanocavity
%           
%           parameter D is given by periods corresponding to 
%           D(1) = d[-n], D(2) = d[-n+1], .. D(n) = d[1] theres no d[0].
% 
function [T_TL_max]=asymmetric_solver_scaled( n, D, lambda, a, h )
   
   if ( nargin < 3 ), error('asymmetric_solver_scaled(n,D,lambda,a,h'), end;
   if ( nargin < 4 ),   a=40,  end;
   if ( nargin < 5 ), h=10.1:2:(lambda/4+5);, end;

   if ( length(D) ~= 2*n )
     error('Too few period lengths for asymmetric periodic solver');
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
   %The following line is needed to obtain e_alpha for a range of depths h   
   %h=100; %to compare with k-domain calculations (single h value), uncomment this line
   for i=1:length(h),
       eps=ones(2*n+1,1)/tan(k*h(i));
       eps(n+1)=-sqrt(-1);
       eps1=diag(eps);
       e_alpha(:,i)=(g-eps1)\rhs;
   end
   transmi_total=1-(abs(1-e_alpha(n+1,:))).^2;

   %Following lines calculate overlap in theta domain (far field, very efficient)
   theta_step=pi/(1000-1);
   theta=(-pi/2):theta_step:(pi/2);
   grid_theta=ones(size(theta));
   grid_theta(1)=1/2;
   grid_theta(length(theta))=1/2;
   h_norm=(e_alpha.'*exp(-sqrt(-1)*k*ds'*sin(theta))*sqrt(k*a/(2*pi))).*(ones(size(h'))*(sin(sin(theta)*k*a/2)./(sin(theta)*k*a/2)));
   z0_step=4*k*a0^2/160; %sign convention: positive z0 is the location of beam waist in free space in front of film.
   z0=((-2*k*a0^2):z0_step:(2*k*a0^2))';
   
   gauss_projection1=((k*a0)^2/pi)^0.25*theta_step*exp(sqrt(-1)*k*z0*cos(theta))*((ones(size(h'))*(exp(-(k*a0*theta).^2/2).*grid_theta)).*h_norm).';
gauss_projection2=(abs(gauss_projection1)).^2;
   
   %The following lines display the results for a range of h values, optimized versus z0
   [t2,in]=max(gauss_projection2);
   optimum_z0=z0(in);
   t=gauss_projection1(in);
   T=t2;
   T_over_L=t2./(transmi_total-t2);
   
   [dummy,idx]=sort(T);%ascending
   [dummy,indx]=sort(T_over_L);%ascending
   idx = idx(end:-1:1); indx = indx(end:-1:1);
   
   T_TL_max.T = T(idx(1));
   T_TL_max.TL = T_over_L(indx(1));
   T_TL_max.n = n;
   T_TL_max.d = D; TL_max.lambda = lambda;
   T_TL_max.h = h(indx(1));
   T_TL_max.z0 = optimum_z0(indx(1));
   T_TL_max.a0 = a0;

end
%
% asymmetric_solver_scaled(2, [300 200 200 250], 324, 40 )
% asymmetric_solver_scaled(4, [100, 200, 300, 400 100, 200, 300, 400], 500, 40 )
%
