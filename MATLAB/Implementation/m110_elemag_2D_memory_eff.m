%% Simulation Parameters
close all; clear; clc

% for getting a better visualization%
% custom color map for values from 0 to 0.5 %
custom_map = [...
    .00 .00 .00 % black
    .20 .00 .00
    .40 .00 .00
    .60 .00 .00
    .80 .00 .00
    1   .00 .00 % red
    1   .20 .00
    1   .40 .00
    1   .60 .00
    1   .80 .00
    1    1  .00 % yellow
    .80  1  .00
    .60  1  .00
    .40  1  .00
    .20  1  .00
    .00  1  .00 % green
    .00 .90 .10
    .00 .80 .20
    .00 .70 .30
    .00 .60 .40
    .00 .50 .50
    .00 .40 .60
    .00 .30 .70
    .00 .20 .80
    .00 .10 .90
    .00 .00  1  % blue
    .05 .00  1
    .10 .00  1
    .15 .00  1
    .20 .00  1
    .25 .00  1
    .30 .00  1
    .35 .00  1
    .40 .00  1
    .45 .00  1
    .50 .00  1
    .55 .00  1
    .60 .00  1
    .65 .00  1
    .70 .00  1
    .75 .00  1
    .80 .00  1
    .85 .00  1
    .90 .00  1
    .95 .00  1
    1   .00  1 % magenta
    1   .05  1
    1   .10  1
    1   .15  1
    1   .20  1
    1   .25  1
    1   .30  1
    1   .35  1
    1   .40  1
    1   .45  1
    1   .50  1
    1   .55  1
    1   .60  1
    1   .65  1
    1   .70  1
    1   .75  1
    1   .80  1
    1   .85  1
    1   .90  1
    1   .95  1
    1    1   1 % white
    ];

%numbero of space and time points
nx = 500;
ny = nx;
nt = 1000;

%constants
eps_0 = 8.854e-12;
eps_1 = eps_0 * 81;
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

%vari�veis auxiliares pra gaussiana
T = 3*(nt*dt/10);
std_dev = 7.005203380146562e-11;

%% Simulate and Video
close all;

% font position %
px1 = 50;          py1 = 400;
px2 = px1 + 1;     py2 = py1;
px3 = px1 + 2;     py3 = py1;

% so the figure won't show up%
figr = figure('Visible', 'On' );

% colormap %
colormap(custom_map)

for n = 0:1:(nt-1)
	t = n*dt;
    
    % font %    
    Ez(px1, py1) = exp(-((t-T).^2)/((std_dev.^2)*.2));
    Ez(px2, py2) = exp(-((t-T).^2)/((std_dev.^2)*.2));
    Ez(px3, py3) = exp(-((t-T).^2)/((std_dev.^2)*.2));
   
    if  (mod(n,2) == 1)
        % plot (not shown)%
        pcolor(Ez)
        shading interp
        colorbar
        caxis([0 .5])
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
            if(j>250)
                eps = eps_1;
            else
                eps = eps_0;
            end
            Ez(i,j) = Ez(i,j) + (dt/eps)*((Hy(i,j+1)-Hy(i,j))/dx - (Hx(i+1,j)-Hx(i,j))/dy);
            
        end
    end
    
    % safety fisrt :D %
    if Ez(px2 +1, py2 +1) > 2
        disp({'Error';['Ez = ', num2str(Ez(px2 +1, py2 +1))];['Iteration: ',num2str(n)]});
        break
    end
end
pause(.1)