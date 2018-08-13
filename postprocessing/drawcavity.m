%%
%% (C) 2007-2009, Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%% This file is part of Plasmonic Nanocavity Solver project 
%% You may use or distribute this file under terms of MIT License
%% 
%% 
%% 03/26/10: Make Matlab 6.0 compatible.
%% 
%% 10/02/08: Draw the geometry of a aperiodic nanocavity 
%%           on current figure. Assumes aperiodic/symmetric
%%           grooves.
%% 
function XYoffset = drawcavity(N,D,h,a,XYoffset)
  true = 1; false = ~true;
  aperiodic = false;

  if ( nargin < 4 )
    return
  end

  ld = length(D);
  if ( ld == 0 )
	error('D length >= 1');
  elseif ( ld == 1)
    aperiodic = false;
  else
    aperiodic = true;
  end

  if ( nargin < 5 )
    XYoffset = [0,0];
  end


  if ( ~aperiodic )
    if ( length(D) == 1)
      D = D*ones(1,N);
    end
  end
  if ( length(h) == 1)
    h = h*ones(1,N);
  end

  

  ht = 2*max(h) ; %total depth.
  X0 = XYoffset(1); Y0 = XYoffset(2) + ht;

  % label this
  text(X0,Y0-ht+5,['N= ',num2str(N),' d(max)=',num2str(max(D)),...
  	      ' d(min)=',num2str(min(D))])
  title(['N= ',num2str(N),' d(max)=',num2str(max(D)),...
  	      ' d(min)=',num2str(min(D))])
  set(gcf,'Color','White')
  % draw through hole trench.
  line([X0 + a/2;X0 + a/2], [Y0 - ht;Y0])
  line([X0 - a/2;X0 - a/2], [Y0 - ht;Y0])


  left = X0 - a/2;  right = X0 + a/2;

  for groove = 1:N
    %% one segment on left with dip and groove.
    line([left,left-D(groove)+a],[Y0,Y0])
    line([left-D(groove)+a,left-D(groove)+a],[Y0,Y0-h(groove)])
    line([left-D(groove)+a,left-D(groove)-a ],[Y0-h(groove) ,Y0-h(groove)])
    line([left-D(groove)-a ,left-D(groove)-a],[Y0-h(groove) ,Y0])
    %% update left
    left = left - D(groove) - a;

    %% one segment on right with dip and groove.
    line([right,right+D(groove)-a],[Y0,Y0])
    line([right+D(groove)-a,right+D(groove)-a],[Y0,Y0-h(groove)])
    line([right+D(groove)-a,right+D(groove)+a] ,[Y0-h(groove),Y0-h(groove)])
    line([right+D(groove)+a ,right+D(groove)+a],[Y0-h(groove) ,Y0])
    %% update right
    right = right + D(groove) + a;
  end


  axis([-16*1e3 +16*1e3  0 Y0])
  XYoffset = [X0,Y0 + max(D) ];
  return

%!test(drawcavity(2,[200 300],[40 40],40,[0 0]))
%!test(drawcavity(4,200,[40 40],40,[0 0]))
%!test(drawcavity(5,150,[40 40],40,[0 0]))
% close all
% clear all
% ofst=drawcavity(2,[200 300],[40 40],40,[0 0]);
% ofst=drawcavity(4,200,50,40,ofst);
% ofst=drawcavity(5,150,[20:10:60],40,ofst);
% axis([-1e3 +1e3  0 ofst(2)])
% pause(2)
% close all
