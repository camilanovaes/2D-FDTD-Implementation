% Clearing variables in memory and Matlab command screen
close all;
clear;
clc;

xdim = 2000;                      % Grid Dimension X
ydim = xdim;                      % Grid Dimension Y
time_tot = 2000;                  % Total of time steps

epsilon0 = 8.854e-12;             % Permittivity of free space (Media 01)
epsilon1 = epsilon0 * 2;          % Permittivity of Media 02
mu0 = pi*4e-7;                    % Permeability of free space
c = 299792458;                    % Speed of electromagnetic wave (m/s)
S = 0.99;                         % Courant stability factor

% Initialization of field matrices
Hx = zeros(xdim,xdim);
Hy = zeros(xdim,xdim);
Ez = zeros(xdim,xdim);

deltax = 1e-4;                                          % Spatial grid step length -> dx: 0.0001 meters
deltay = deltax;                                        % Spatial grid step length -> dy: 0.0001 meters
deltat = (S/(c*sqrt(1/(deltax.^2)+1/(deltay.^2))));     % Temporal grid step obtained using Courant condition

% Gaussian font parameters
T = (time_tot*deltat/10);

% Initialization of horn antenna
Ax = 750;
Ay = 500;
BoxLeftSide = zeros(5,1);
BoxTop = zeros(1,6);
BoxBottom = zeros(1,8);
TopStep = ones(9,10);
TopDiag = zeros(10,1);
BottomDiag = ones(5,5);

for i = 1:1:10
  if (i > 1 && i ~= 10)
    TopStep(i,i) = 0;
    TopStep(i-1,i) = 0;
  elseif(i == 10)
    TopStep(i-1,i) = 0;
  else
    TopStep(i,i) = 0;
  end
end
BottomStep = TopStep;

for j = 1:1:5
  BottomDiag(j,j) = 0;
end

% Font position
xsource1 = Ax+1;      
ysource1 = Ay+5;
xsource2 = xsource1 + 1;     
ysource2 = ysource1;
xsource3 = xsource1 + 2;
ysource3 = ysource1;

figr = figure('Visible', 'On' );
colormap('hsv');

for n = 0:1:time_tot
  t = n*deltat;
    
  % Font    
  Ez(xsource1, ysource1) = exp(-((t-(3*T))^2)/((T^2)))* sin(2*pi*30e9*t);
  Ez(xsource2, ysource2) = exp(-((t-(3*T))^2)/((T^2)))* sin(2*pi*30e9*t);
  Ez(xsource3, ysource3) = exp(-((t-(3*T))^2)/((T^2)))* sin(2*pi*30e9*t);
   
  if (mod(n,2) == 1)
    pcolor(-log(abs(Ez)+1e-30))
    shading interp
    colorbar
    title({['Nt = ',num2str(n)],['Time: ',num2str(t),' sec.']});
    hold on;
    plot ([ysource2, ydim-xsource2+ysource2],[xsource2, ydim],'k','linewidth',1);
    hold off
    pause(.005)
  end
  pause(.005)
    
  for i = 2:1:xdim-1
    for j = 2:1:ydim-1
      Hx(i,j) = Hx(i,j) - (deltat/mu0)*(Ez(i,j)-Ez(i-1,j))/deltay;
      Hy(i,j) = Hy(i,j) + (deltat/mu0)*(Ez(i,j)-Ez(i,j-1))/deltax;
    end
  end
    
  for i = 1:1:xdim-1
    for j = 1:1:ydim-1
      if(i>xdim/2)
        eps = epsilon1;                 % Medium 02
      elseif (i == xdim/2)
        eps = (epsilon0 + epsilon1)/2;  % Middle Interface
      else
        eps = epsilon0;                 % Medium 01
      end 
      Ez(i,j) = Ez(i,j) + (deltat/eps)*((Hy(i,j+1)-Hy(i,j))/deltax - (Hx(i+1,j)-Hx(i,j))/deltay);        
    end
  end
    
  % Boundary Conditions for horn antenna
  Ez(Ax:1:(Ax)+4,Ay) = BoxLeftSide;
  Ez((Ax)-1,(Ay):(Ay)+7) = BoxBottom;
  Ez((Ax)+5,(Ay):(Ay)+5) = BoxTop;
  Ez((Ax)+6:(Ax)+14,(Ay)+5:(Ay)+14) = Ez((Ax)+6:(Ax)+14,(Ay)+5:(Ay)+14).*TopStep;
  Ez((Ax):(Ax)+8,(Ay)+7:(Ay)+16) = Ez((Ax):(Ax)+8,(Ay)+7:(Ay)+16).*BottomStep;
  Ez((Ax)+15:(Ax)+24,(Ay)+14) = TopDiag;
  Ez((Ax)+7:-1:(Ax)+3,(Ay)+17:(Ay)+21) = Ez((Ax)+7:-1:(Ax)+3,(Ay)+17:(Ay)+21).*BottomDiag;
    
end
pause(.1)
