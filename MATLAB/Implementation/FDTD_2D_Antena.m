%% Simulation Parameters
close all;
clear;
clc

%numbero of space and time points
xdim = 1500;
ydim = xdim;
time_tot = 1000;

%constants
epsilon0 = 8.854e-12;
epsilon1 = epsilon0 * 2;
mu0 = pi*4e-7;
c = 299792458; 	% m/s
S = 0.99; %Courant stability factor

%fields
Hx = zeros(xdim,xdim);
Hy = zeros(xdim,xdim);
Ez = zeros(xdim,xdim);

%diferenctial elements
deltax = 1e-3;              % 0.001 meters
deltay = deltax;
deltat = (S/(c*sqrt(1/(deltax.^2)+1/(deltay.^2))));     %formula de estabilidade

%variáveis auxiliares pra gaussiana
T = 3*(time_tot*deltat/10);
std_dev = 7.005203380146562e-11;

%Inicializando antena
Ax = 220;
Ay = 100;
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

% font position %
xsource1 = Ax+1;      
ysource1 = Ay+5;
xsource2 = xsource1 + 1;     
ysource2 = ysource1;
xsource3 = xsource1 + 2;
ysource3 = ysource1;

% settings to save video %
orig_file = 'D:\Desktop\eletromag\';
cd 'D:\Desktop\eletromag\';
vid_obj = VideoWriter('FDTD_2D_Antena.avi');
vid_obj.FrameRate = 30;
cd (orig_file);
%% Simulate and Video
close all;
% so the figure won't show up%
figr = figure('Visible', 'Off' );

% wait bar %
wb = waitbar(0,'Please wait...','Name','Calculating and Generating video');

% open video object%
open(vid_obj);

%colormap('jet');
%colormap('colorcube');
%colormap('prism');
colormap('hsv');

for n = 0:1:time_tot
	
    t = n*deltat;
    
    % font %    
    Ez(xsource1, ysource1) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
    Ez(xsource2, ysource2) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
    Ez(xsource3, ysource3) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
   
    if  (mod(n,2) == 1)
        % plot (not shown)%
        pcolor(-log(abs(Ez)+1e-30))
        shading interp
        colorbar
        title({['Nt = ',num2str(n)],['Time: ',num2str(t),' sec.']});
		hold on;
        plot ([ysource2, ydim-xsource2+ysource2],[xsource2, ydim],'k','linewidth',1);
		hold off
        pause(.005)
        % writing video from the figure %
        img = getframe(figr);
        writeVideo(vid_obj, img);
    end
    pause(.005)
     % wait bar %
    waitbar(n/time_tot, wb, {'Please wait...';['Nt = ',num2str(n),' /',num2str(time_tot)]});
    
    for i = 2:1:xdim-1
        for j = 2:1:ydim-1
          Hx(i,j) = Hx(i,j) - (deltat/mu0)*(Ez(i,j)-Ez(i-1,j))/deltay;
          Hy(i,j) = Hy(i,j) + (deltat/mu0)*(Ez(i,j)-Ez(i,j-1))/deltax;
        end
    end
    
    for i = 1:1:xdim-1
        for j = 1:1:ydim-1
            if(i>xdim/2)
                eps = epsilon1;
            elseif (i == xdim/2)
                eps = (epsilon0 + epsilon1)/2;
            else
                eps = epsilon0;
            end 
                Ez(i,j) = Ez(i,j) + (deltat/eps)*((Hy(i,j+1)-Hy(i,j))/deltax - (Hx(i+1,j)-Hx(i,j))/deltay);        
        end
    end
    
    % colocar a antena
    Ez(Ax:1:(Ax)+4,Ay) = BoxLeftSide;
    Ez((Ax)-1,(Ay):(Ay)+7) = BoxBottom;
    Ez((Ax)+5,(Ay):(Ay)+5) = BoxTop;
    Ez((Ax)+6:(Ax)+14,(Ay)+5:(Ay)+14) = Ez((Ax)+6:(Ax)+14,(Ay)+5:(Ay)+14).*TopStep;
    Ez((Ax):(Ax)+8,(Ay)+7:(Ay)+16) = Ez((Ax):(Ax)+8,(Ay)+7:(Ay)+16).*BottomStep;
    Ez((Ax)+15:(Ax)+24,(Ay)+14) = TopDiag;
    Ez((Ax)+7:-1:(Ax)+3,(Ay)+17:(Ay)+21) = Ez((Ax)+7:-1:(Ax)+3,(Ay)+17:(Ay)+21).*BottomDiag;
end
waitbar(1, wb, {'Please wait...';'Saving...'});

% close and finish video%
close(vid_obj);
pause(.1)
% close wait bar%
delete(wb)
