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

%constants
eps_0 = 8.854e-12;
mi_0 = pi*4e-7;
C_0 = 299792458; 	% m/s

%numbero of space and time points
nx = 500;
ny = nx;
nt = 1000;

%fields
Ex = zeros(nx,nx);
Ey = zeros(nx,nx);
Hz = zeros(nx,nx);

%diferenctial elements
dx = 1e-3;              % 0.001 meters
dy = dx;
dt = (0.99/(C_0*sqrt(1/(dx.^2)+1/(dy.^2))));     %formula de estabilidade

%variáveis auxiliares pra gaussiana
T = 3*(nt*dt/10);
std_dev = 7.005203380146562e-11;

% settings to save video %
orig_file = 'D:\Desktop\eletromag\';
cd 'D:\Desktop\eletromag\';
vid_obj = VideoWriter('m110_elemag_2D_memory_eff.avi');
vid_obj.FrameRate = 30;
cd (orig_file);

%% Simulate and Video
close all;

% font position %
px1 = 50;          py1 = 50;
px2 = px1 + 1;     py2 = py1;
px3 = px1 + 2;     py3 = py1;

% so the figure won't show up%
figr = figure('Visible', 'Off' );

% colormap %
colormap(custom_map)

% wait bar %
wb = waitbar(0,'Please wait...','Name','Calculating and Generating video');

% open video object%
open(vid_obj);

for n = 0:1:(nt-1)
	t = n*dt;
    
    % font %    
    Hz(px1, py1) = exp(-((t-T).^2)/((std_dev.^2)*.2));
    Hz(px2, py2) = exp(-((t-T).^2)/((std_dev.^2)*.2));
    Hz(px3, py3) = exp(-((t-T).^2)/((std_dev.^2)*.2));
   
    if  (mod(n,2) == 1)
        % plot (not shown)%
        pcolor(Hz)
        shading interp
        colorbar
        caxis([0 .5])
        title({['Nt = ',num2str(n)],['Time: ',num2str(t),' sec.']});
        
        pause(.005)
        
        % writing video from the figure %
        img = getframe(figr);
        writeVideo(vid_obj, img);
    end
    pause(.005)
    
    % wait bar %
    waitbar(n/(nt-1), wb, {'Please wait...';['Nt = ',num2str(n),' /',num2str(nt-1)]});
    
    % antenna %
    %      (py-2) (py1)
    %          |   |
    %          v   v      (py+4)
    % px1 -1 ->N N N N N   |
    % px1    ->N   O   N   v  (py+6)
    % px1 +1 ->N   o   N N     |
    % px1 +2 ->N   o     N N   v  (py+8)
    %          N N         N N     |
    % px1 +4 ->  N N         N N   v  (py+10)
    %              N N         N N     |
    %     px1 +6 ->  N N         N N   v
    %                  N N         N N        (py+14)
    %       px1 +8 ->    N N         N N       |
    %                      N N         N N     v
    %           px1 +10 ->   N N         N N N N <- px1 +10
    %                          N N
    %             px1 +12 ->     N N
    %                              N
    %                 px1 +14 ->   N
    %                 px1 +15 ->   N
    
    
    % maxwell and FDTD stuff %
    for i = 2:1:nx-1
        for j = 2:1:ny-1
            
            % if to limit the antenna %
            if (    (i ~= (px1 -1) || j ~= (py1 -2)) && ...
                    (i ~= (px1 -1) || j ~= (py1 -1)) && ...
                    (i ~= (px1 -1) || j ~= (py1   )) && ...
                    (i ~= (px1 -1) || j ~= (py1 +1)) && ...
                    (i ~= (px1 -1) || j ~= (py1 +2)) && ...
                    ...
                    (i ~= (px1   ) || j ~= (py1 -2)) && ...
                    (i ~= (px1   ) || j ~= (py1 +2)) && ...
                    ...
                    (i ~= (px1 +1) || j ~= (py1 -2)) && ...
                    (i ~= (px1 +1) || j ~= (py1 +2)) && ...
                    (i ~= (px1 +1) || j ~= (py1 +3)) && ...
                    ...
                    (i ~= (px1 +2) || j ~= (py1 -2)) && ...
                    (i ~= (px1 +2) || j ~= (py1 +3)) && ...
                    (i ~= (px1 +2) || j ~= (py1 +4)) && ...
                    ...
                    (i ~= (px1 +3) || j ~= (py1 -2)) && ...
                    (i ~= (px1 +3) || j ~= (py1 -1)) && ...
                    (i ~= (px1 +3) || j ~= (py1 +4)) && ...
                    (i ~= (px1 +3) || j ~= (py1 +5)) && ...
                    ...
                    (i ~= (px1 +4) || j ~= (py1 -1)) && ...
                    (i ~= (px1 +4) || j ~= (py1   )) && ...
                    (i ~= (px1 +4) || j ~= (py1 +5)) && ...
                    (i ~= (px1 +4) || j ~= (py1 +6)) && ...
                    ...
                    (i ~= (px1 +5) || j ~= (py1   )) && ...
                    (i ~= (px1 +5) || j ~= (py1 +1)) && ...
                    (i ~= (px1 +5) || j ~= (py1 +6)) && ...
                    (i ~= (px1 +5) || j ~= (py1 +7)) && ...
                    ...
                    (i ~= (px1 +6) || j ~= (py1 +1)) && ...
                    (i ~= (px1 +6) || j ~= (py1 +2)) && ...
                    (i ~= (px1 +6) || j ~= (py1 +7)) && ...
                    (i ~= (px1 +6) || j ~= (py1 +8)) && ...
                    ...
                    (i ~= (px1 +7) || j ~= (py1 +2)) && ...
                    (i ~= (px1 +7) || j ~= (py1 +3)) && ...
                    (i ~= (px1 +7) || j ~= (py1 +8)) && ...
                    (i ~= (px1 +7) || j ~= (py1 +9)) && ...
                    ...
                    (i ~= (px1 +8) || j ~= (py1 +3)) && ...
                    (i ~= (px1 +8) || j ~= (py1 +4)) && ...
                    (i ~= (px1 +8) || j ~= (py1 +9)) && ...
                    (i ~= (px1 +8) || j ~= (py1 +10)) && ...
                    ...
                    (i ~= (px1 +9) || j ~= (py1 +4)) && ...
                    (i ~= (px1 +9) || j ~= (py1 +5)) && ...
                    (i ~= (px1 +9) || j ~= (py1 +10)) && ...
                    (i ~= (px1 +9) || j ~= (py1 +11)) && ...
                    ...
                    (i ~= (px1 +10) || j ~= (py1 +5)) && ...
                    (i ~= (px1 +10) || j ~= (py1 +6)) && ...
                    (i ~= (px1 +10) || j ~= (py1 +11)) && ...
                    (i ~= (px1 +10) || j ~= (py1 +12)) && ...
                    (i ~= (px1 +10) || j ~= (py1 +13)) && ...
                    (i ~= (px1 +10) || j ~= (py1 +14)) && ...
                    ...
                    (i ~= (px1 +11) || j ~= (py1 +6)) && ...
                    (i ~= (px1 +11) || j ~= (py1 +7)) && ...
                    ...
                    (i ~= (px1 +12) || j ~= (py1 +7)) && ...
                    (i ~= (px1 +12) || j ~= (py1 +8)) && ...
                    ...
                    (i ~= (px1 +13) || j ~= (py1 +8)) && ...
                    ...
                    (i ~= (px1 +14) || j ~= (py1 +8)) && ...
                    ...
                    (i ~= (px1 +15) || j ~= (py1 +8)) )
                
                Ex(i,j) = Ex(i,j) + (dt/eps_0)*(Hz(i,j)-Hz(i-1,j))/dy;
                
                Ey(i,j) = Ey(i,j) - (dt/eps_0)*(Hz(i,j)-Hz(i,j-1))/dx;
            end
        end
    end
    
    % more maxwell and FDTD stuff %
    for i = 1:1:nx-1
        for j = 1:1:ny-1
            Hz(i,j) = Hz(i,j) - (dt/mi_0)*((Ey(i,j+1)-Ey(i,j))/dx - (Ex(i+1,j)-Ex(i,j))/dy);
        end
    end
    
    % safety fisrt :D %
    if Hz(px2 +1, py2 +1) > 2
        disp({'Error';['Hz = ', num2str(Hz(px2 +1, py2 +1))];['Iteration: ',num2str(n)]});
        break
    end
end

waitbar(1, wb, {'Please wait...';'Saving...'});

% close and finish video%
close(vid_obj);

pause(.1)

% close wait bar%
delete(wb)






