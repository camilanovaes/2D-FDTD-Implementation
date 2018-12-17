%Clearing variables in memory and Matlab command screen

clear all;
clc;

% Grid Dimension in x (xdim) and y (ydim) directions
xdim=200;
ydim=200;

%Total no of time steps
time_tot=350;

%Position of the source (center of the domain)
xsource=100;
ysource=100;

%Courant stability factor
S = 0.99;

% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0 = 8.854*1e-12;
mu0 = 4*pi*1e-7;
c = 2.99e+8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IMPORTANTISSIMO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------%
% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta = 1e-6;
% Temporal grid step obtained using Courant condition
deltat = S*delta/c;
%---------------------------------------------------------------------------%

% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

% Initialization of permittivity and permeability matrices
epsilon = epsilon0;
mu = mu0;

%Choice of nature of source
gaussian=1;
sine=0;

% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
impulse=0;

% Update loop begins
for n=1:1:time_tot
    
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource,ysource)=1;
    end
    
    %Vector update instead of for-loop for Hy and Hx fields
    for l = 1: xdim 
        for m = 1:ydim-1
            Hx(m,l) = Hx(m,l) - 0.5*(Ez(m+1,l)-Ez(m,l));
        end
    end
    
    for m1 = 1: ydim
        for l1 = 1:xdim-1
            Hy(m1,l1) = Hy(m1,l1) + 0.5*(Ez(m1,l1+1)-Ez(m1,l1));
        end
    end
    
    for m2 = 2: ydim
        for l2 = 2:xdim
            Ez(m2,l2) = Ez(m2,l2) + 0.5*(Hy(m2,l2) - Hy(m2,l2-1)) - 0.5*(Hx(m2,l2)-Hx(m2-1,l2));
        end
    end
    
    % Perfect Electric Conductor boundary condition
    Ez(1:xdim,1)=0;
    Ez(1:xdim,ydim)=0;
    Ez(1,1:ydim)=0;
    Ez(xdim,1:ydim)=0;
    
    % Source conditions
    if impulse==0
        % If unit-time step
        if gaussian==0 && sine==0
            Ez(xsource,ysource)=1;
        end
        %if sine
        if sine==1
            tstart=1;
            N_lambda=c/(frequency*delta);
            Ez(xsource,ysource)=sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        %if gaussian
        if gaussian==1
            if n<=42
                Ez(xsource,ysource)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
                Ez(xsource,ysource+1)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
                Ez(xsource,ysource+2)=(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
            else
                Ez(xsource,ysource)=0;
            end
        end
    else
        %if impulse
        Ez(xsource,ysource)=0;
    end
    
    %Movie type colour scaled image plot of Ez
    imagesc(delta*(1:1:xdim)*1e+6,(1e+6*delta*(1:1:ydim))',Ez',[-1,1]);colorbar;
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with PEC boundary and at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
