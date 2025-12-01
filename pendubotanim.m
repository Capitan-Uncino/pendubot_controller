function [sys,x0]=pendubotanim(t,x,u,flag,ts);

%||=================================================||
%||                                                 ||
%||   INMA 2671                                     ||
%||                                                 ||
%||   Pr J. Hendrickx, Eng. A. Taylor & F. Wielant  ||
%||                                                 ||
%||=================================================||

%pendubotanim S-function for animating the motion of the Pendubot.

%   Ned Gulley, 6-21-93
%   Copyright (c) 1990-1998 by The MathWorks, Inc. All Rights Reserved.
%   $Revision: 5.6 $
%   Editted by Dan Block

global PendAnim2

if flag==2,
  
  shh = get(0,'ShowHiddenHandles');
  set(0,'ShowHiddenHandles','on');

  if any(get(0,'Children')==PendAnim2),
    if strcmp(get(PendAnim2,'Name'),'Pendubot Animation'),
      set(0,'currentfigure',PendAnim2);
      hndl=get(gca,'UserData');
      b=1; c=1;
      xSFP=0; ySFP=0;
      xBFP=xSFP+b*sin(u(1)); yBFP=ySFP-0.85*b*cos(u(1));
      xMC=xBFP+c*sin(u(2)); yMC=yBFP-c*cos(u(2));
      x=[xSFP xBFP NaN xBFP xMC];
      y=[ySFP yBFP NaN yBFP yMC];
      %set(hndl,'XData',x,'YData',y);
      set(hndl,'XData',x,'YData',y, ...
          'EraseMode','normal', ...
           'LineWidth',7, ...
           'Marker','o', ...
           'MarkerFaceColor' , 'k' , ...
           'MarkerEdgeColor' , 'r' , ...
           'MarkerSize',9);

      drawnow;
    end
  end
    
  set(0,'ShowHiddenHandles',shh);
  
  sys=[];

elseif flag==0,
  % Initialize the figure for use with this simulation
  [PendAnim2 PendAnim2Axes] = animinit('Pendubot Animation');	  
  %animinit('Pendubot Animation');
  %[flag,PendAnim2] = figflag('Pendubot Animation');
   if ishghandle(PendAnim2)
      % bring figure to foreground
      figure(PendAnim2);
  end
   
 
  axis(PendAnim2Axes,[-3 3 -2 2]);
  hold(PendAnim2Axes, 'on'); 
  %axis([-2 2 -3 2]);
  %hold on;

     % add the body
  legh = 0.6;
  bodyw = 2;
  legw = 2;
  bodyh = 0;
  mount = 0.6;
  set(0,'currentfigure',PendAnim2);   
  plot(PendAnim2Axes, [-legw/2,-legw/2,-legw/2,...
      -bodyw/2,-bodyw/2,bodyw/2,bodyw/2,...
      legw/2,legw/2,legw/2,-legw/2],...
      -.1+[0,-legh,0,0,bodyh,bodyh,0,0,-legh,0,0]...
      -(bodyh+mount/2),'linestyle','-','color',[.4 .4 .4],'LineWidth',7);
  plot(PendAnim2Axes, [-mount/2,-mount/2,mount/2,mount/2],...
      [-mount/2,mount/2,mount/2,-mount/2],...
      'linestyle','-','color',[0 0 0],'LineWidth',2);
  
  % Set up the geometry for the problem
  % SFP=Space Fixed Pivot
  % BFP=Body Fixed Pivot
  b=1; c=1;
  xSFP=0; ySFP=0;
  xBFP=xSFP; yBFP=ySFP-0.85*b;
  xMC=xBFP; yMC=yBFP-c;
  x=[xSFP xBFP NaN xBFP xMC];
  y=[ySFP yBFP NaN yBFP yMC];
  hndl=plot(PendAnim2Axes,x,y, ...
           'EraseMode','background', ...
           'LineWidth',8, ...
           'Marker','.', ...
           'MarkerSize',18);
  %set(gca,'DataAspectRatio',[1 1 1]);
  %set(gca,'UserData',hndl);
  set(PendAnim2Axes,'DataAspectRatio',[1 1 1]);
  set(PendAnim2Axes,'UserData',hndl);

  sys=[0 0 0 2 0 0];
  x0=[];
  str = [];
  ts  = [-1, 0];
  % specify that the simState for this s-function is same as the default
  simStateCompliance = 'DefaultSimState';

end

