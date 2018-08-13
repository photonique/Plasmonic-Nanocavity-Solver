%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%
% 04/27/09: If R -> imaginary values we invoke discard_flag=true; for
%           that geometry.
% 04/16/09: Include the R -> 0.99max in Ecav/Ein calculation.
% 04/15/09: Limit the R to 0.99 max; suggested by Dr.V.
% 04/09/09: Add R parameter as well.
% 04/08/09: Add Ttotal and cavity enhancement factors for the geometry
%           at the given wavelength:
%           Ecav/Ein = t/(1 - r) = sqrt(Tgauss)/(1 - sqrt(1 - Ttotal));
% 
% 10/18/08: Certain choice of parameters lend to singular matrices
%               and we necessarily need to change the solver for
%                inverting the matrix. Solve using LU decomposition.
%                Add a discard flag for the data.
% 
% 10/08/08: solver extended to span the variable groove depth.
%           user supplies hrange_cartprod as cartesian prod depths to optimize
% 12/19/07: solver for aperiodic case.
% 12/15/07: convert solver for use with a GA.
% 
function [T_TL_max,discard_flag]=aperiodic_depthsolver_scaled( n, D, lambda, a, hrange_cartprod )
   
   if ( nargin < 3 ), error('aperiodic_solver_scaled_nanocavity'), end;
   if ( nargin < 4 ),   a=40,  end;
   if ( nargin < 5 ), h=10.1:2:(lambda/4+5);, end;

if ( 0 )
   n
   lambda
   D
   hrange_cartprod 
end

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
	T_TL_max={};
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
	 T_TL_max={};
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
   transmi_total=1-(abs(1-e_alpha(n+1,:))).^2;
   
   %Following lines calculate overlap in theta domain (far field, very efficient)
   theta_step=pi/(1000-1); 
   theta=(-pi/2):theta_step:(pi/2); 
   grid_theta=ones(size(theta)); 
   grid_theta(1)=1/2; 
   grid_theta(length(theta))=1/2;
   
   h_norm=(e_alpha.'*exp(-sqrt(-1)*k*ds'*sin(theta))*sqrt(k*a/(2*pi))).*(ones([rows(hrange_cartprod),1])*(sin(sin(theta)*k*a/2)./(sin(theta)*k*a/2))); 
   z0_step=4*k*a0^2/160; %sign convention: positive z0 is the location of beam waist in free space in front of film. 
   z0=((-2*k*a0^2):z0_step:(2*k*a0^2))';
   
   gauss_projection1=((k*a0)^2/pi)^0.25*theta_step*exp(sqrt(-1)*k*z0*cos(theta))*(((ones([rows(hrange_cartprod),1]))*(exp(-(k*a0*theta).^2/2).*grid_theta)).*h_norm).';
   gauss_projection2=(abs(gauss_projection1)).^2; 
    
   %The following lines display the results for a range of h values, optimized versus z0 
   [t2,in]=max(gauss_projection2); 
   optimum_z0=z0(in); 
   t=gauss_projection1(in); 
   T=t2; 
   T_over_L=t2./(transmi_total-t2); 
    
   [dummy,idx]=sort(T,'descend'); 
   [dummy,indx]=sort(T_over_L,'descend');

   Tmax = T(idx(1));  TLmax = T_over_L(indx(1));   
   Ttotal = Tmax*( 1 + 1/TLmax);
   R = 0.99 - Ttotal;
   
   if ( R < 0 )
	T_TL_max={};
	discard_flag=true;
	n
	D
	lambda
	hrange_cartprod
	printf('R < 0. Failed!\n')
	return
   end
   
   Ecav_Ein = sqrt( Tmax )/( 1 - sqrt( R ) );
   
   T_TL_max.T = Tmax;
   T_TL_max.TL = TLmax;
   T_TL_max.Ttotal = Ttotal;
   T_TL_max.R = R;
   T_TL_max.Ecav_Ein = Ecav_Ein;
   T_TL_max.n = n;
   T_TL_max.d = D; TL_max.lambda = lambda;
   T_TL_max.h = hrange_cartprod(indx(1),:);
   T_TL_max.z0 = optimum_z0(indx(1));
   T_TL_max.a0 = a0;
   ## T_TL_max.lambda = lambda; ##? required?

end
%
% aperiodic_depthsolver_scaled(2, [300 300], 324, [40 50] )
% aperiodic_depthsolver_scaled(4, [100, 200, 300, 400], 500,[45 40 45 40] )
%
