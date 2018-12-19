%% Simulation Parameters
close all; clear; clc

%numbero of space and time points
nx = 100;
ny = nx;
nt = 1000;

%constants
eps_0 = 8.854e-12;
eps_1 = eps_0 * 5;
mi_0 = pi*4e-7;
C_0 = 299792458; 	% m/s

%fields
Hx = zeros(nx,nx);
Hy = zeros(nx,nx);
Ez = zeros(nx,nx);

%diferenctial elements
dx = 1e-3;              % 0.001 meters
dy = dx;
dt = (0.99/(C_0*sqrt(1/(dx.^2)+1/(dy.^2))));     %formula de estabilidade

%variáveis auxiliares pra gaussiana
T = 3*(nt*dt/10);
std_dev = 7.005203380146562e-11;

%% Simulate and Video
close all;

% font position %
px1 = nx-94;      py1 = ny-90;
px2 = px1 + 1;     py2 = py1;
px3 = px1 + 2;     py3 = py1;

%Inicializando antena
BoxLeftSide = zeros(5,1);
BoxTop = zeros(1,6);
BoxBottom = zeros(1,8);
TopStep = ones(9,10);
TopDiag = ones(11,4);

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

% so the figure won't show up%
figr = figure('Visible', 'On' );

colormap('jet');

for n = 0:1:(nt-1)
	t = n*dt;
    
    % font %    
    Ez(px1, py1) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
    Ez(px2, py2) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
    Ez(px3, py3) = exp(-((t-T).^2)/((std_dev.^2)*.2))* sin(2*pi*30e9*t);
   
    if  (mod(n,2) == 1)
        % plot (not shown)%
        pcolor(-log(abs(Ez)+1e-30))
        shading interp
        colorbar
        %caxis([0 .5])
        title({['Nt = ',num2str(n)],['Time: ',num2str(t),' sec.']});
        
        pause(.005)
    end
    pause(.005)
    
    for i = 2:1:nx-1
        for j = 2:1:ny-1
                
                Hx(i,j) = Hx(i,j) - (dt/mi_0)*(Ez(i,j)-Ez(i-1,j))/dy;
                Hy(i,j) = Hy(i,j) + (dt/mi_0)*(Ez(i,j)-Ez(i,j-1))/dx;
        end
    end
    
    for i = 1:1:nx-1
        for j = 1:1:ny-1
            if(i>50)
                eps = eps_1;
            elseif (i == 50)
                eps = (eps_0 + eps_1)/2;
            else
                eps = eps_0;
            end 
                Ez(i,j) = Ez(i,j) + (dt/eps)*((Hy(i,j+1)-Hy(i,j))/dx - (Hx(i+1,j)-Hx(i,j))/dy);
        end
    end
    
    % colocar a antena
    Ez(nx-95:1:(nx-95)+4,ny-95) = BoxLeftSide;
    Ez((nx-95)-1,(ny-95):(ny-95)+7) = BoxBottom;
    Ez((nx-95)+5,(ny-95):(ny-95)+5) = BoxTop;
    Ez((nx-95)+6:(nx-95)+14,(ny-95)+5:(ny-95)+14) = Ez((nx-95)+6:(nx-95)+14,(ny-95)+5:(ny-95)+14).*TopStep;
    Ez((nx-95):(nx-95)+8,(ny-95)+7:(ny-95)+16) = Ez((nx-95):(nx-95)+8,(ny-95)+7:(ny-95)+16).*BottomStep;
    
    % safety fisrt :D %
    if Ez(px2 +1, py2 +1) > 2
        disp({'Error';['Ez = ', num2str(Ez(px2 +1, py2 +1))];['Iteration: ',num2str(n)]});
        break
    end
end
pause(.1)