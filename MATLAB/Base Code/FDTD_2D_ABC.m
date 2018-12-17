%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "2D FDTD solution for Mur's Absorbing Boundary Condition (ABC)"
% 
% Objective of the program is to solve for the Maxwell's equation for a TM 
% wave containing the xy-plane polarized magnetic field having components Hy
% and Hx and z-polarized electric field Ez. The fields are updated at every 
% timestep, in a space, where all physical parameters of free space are not
% normalized to 1 but are given real and known values. The update is done 
% using standard update equations obtained from the difference form of Maxwell's 
% curl equations with very low electric and magnetic conductivities of 4x10^(-4) 
% units incorporated in them. The field points are defined in a grid described 
% by Yee's algorithm. The H fields are defined at every half coordinate of 
% spacesteps. More precisely, the Hx part is defined at every half y-coordinate 
% and full x-coordinate and the Hy part is defined at every half x-coordinate 
% and full y-coordinate and E fields i.e the Ez part is defined at every full 
% x and full y-coordinate points.Also here, the space-step length is taken 
% as 1 micron instead of 1 unit in unitless domain assumed in previous programs. 
% Also, the time update is done using Leapfrog time-stepping. Here, H-fields
% i.e. Hx and Hy are updated every half time-step and E fields i.e Ez are 
% updated every full time-step. This is shown by two alternating vector updates 
% spanning only a part of spatial grid where the wave, starting from source, 
% has reached at that particular time instant avoiding field updates at all 
% points in the grid which is unnecessary at that time instant. These spatial 
% updates are inside the main for-loop for time update, spanning the entire 
% time grid. Also, here, the matrices used as multiplication factors for update 
% equations are initialized before the loop starts to avoid repeated calculation 
% of the same in every loop iteration, a minor attempt at optimization. The 
% boundary condition here is Mur's Absorbing Boundary Condition (ABC) where 
% the fields at the grid points have electric field values formulated using 
% Engquist Majda one way wave equations [1] where the boundaries give a sense 
% of absorbing the total field incident on them and reflecting none back to 
% the domain.
%
% A source of electric field is defined at the center of the spatial domain, 
% which is a hard source, in that it does not change its value due to 
% interference from external fields i.e in other words, the source is a 
% perfect electric conductor. The form of the source can be varied using 
% the variables sine, gaussian and impulse. The source is available in four 
% standard forms- Unit-time step, Impulse, Gausian and Sinusoidal forms. 
% The color scaled plot of Ez field over the entire spatial domain is shown 
% at every time step. The simulation can be ended by closing this plot window 
% or by waiting till all the time step updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
S=1/(2^0.5);
% Parameters of free space (permittivity and permeability and speed of
% light) are all not 1 and are given real values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;
% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;
% Initialization of field matrices
Ez=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);
% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);
% Initializing electric and magnetic conductivity matrices
sigma=4e-4*ones(xdim,ydim);
sigma_star=4e-4*ones(xdim,ydim);
%Choice of nature of source
gaussian=0;
sine=0;
% The user can give a frequency of his choice for sinusoidal (if sine=1 above) waves in Hz 
frequency=1.5e+13;
impulse=0;
%Choose any one as 1 and rest as 0. Default (when all are 0): Unit time step
%Multiplication factor matrices for H matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed
A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star)); 
B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
                          
%Multiplication factor matrices for E matrix update to avoid being calculated many times 
%in the time update loop so as to increase computation speed                          
C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma)); 
D=(deltat/delta)./(epsilon+0.5*deltat*sigma);  
%Mur's absorbing boundary condition parameters
p0=1;
p2=-0.5;
%Co-efficients of present and previous (regarding time-step) boundary Ez values 
%in boundary update equation for forward/up and backward/down boundaries in 
%the domain (as given in [1])
c0=(c/(2*S))*(1-(p0/S));
c1=-(c/(2*S))*(1+(p0/S));
c2=(c/(S^2))*(p0+(p2*S*S));
c3=-(p2*c)/2;
c0efffor=-(c0/c1);
c2efffor=-(c2/c1);
c3efffor=-(c3/c1);
c0=(c/(2*S))*(1+(p0/S));
c1=-(c/(2*S))*(1-(p0/S));
c2=-(c/(S^2))*(p0+(p2*S*S));
c3=(p2*c)/2;
c1effrev=-(c1/c0);
c2effrev=-(c2/c0);
c3effrev=-(c3/c0);
%Storage vectors for Ez boundary and boundary-1 values of previous and its
%previous timesteps
prev_xfor=zeros(1,ydim);
prev_x_minus_1for=zeros(1,ydim);
prev_yfor=zeros(xdim,1);
prev_y_minus_1for=zeros(xdim,1);
prev_xrev=zeros(1,ydim);
prev_x_minus_1rev=zeros(1,ydim);
prev_yrev=zeros(xdim,1);
prev_y_minus_1rev=zeros(xdim,1);
% Update loop begins
for n=1:1:time_tot
    
    %if source is impulse or unit-time step 
    if gaussian==0 && sine==0 && n==1
        Ez(xsource,ysource)=1;
    end
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % vector where the wave has reached to avoid unnecessary updates.
    if n<xsource-3
        n1=xsource-n-1;
    else
        n1=2;
    end
    if n<xdim-2-xsource
        n2=xsource+n;
    else
        n2=xdim-2;
    end
    if n<ysource-3
        n11=ysource-n-1;
    else
        n11=2;
    end
    if n<ydim-2-ysource
        n22=ysource+n;
    else
        n22=ydim-2;
    end
    
    %Vector update instead of for-loop for Hy and Hx fields
    Hx(n1:n2-1,n11:n22-1)=A(n1:n2-1,n11:n22-1).*Hx(n1:n2-1,n11:n22-1)-B(n1:n2-1,n11:n22-1).*(Ez(n1:n2-1,n11+1:n22)-Ez(n11:n2-1,n11:n22-1));
    Hy(n1:n2-1,n11:n22-1)=A(n1:n2-1,n11:n22-1).*Hy(n1:n2-1,n11:n22-1)+B(n1:n2-1,n11:n22-1).*(Ez(n1+1:n2,n11:n22-1)-Ez(n11:n2-1,n11:n22-1));
    
    %Vector update instead of for-loop for Ez field
    Ez(n1+1:n2-1,n11+1:n22-1)=C(n1+1:n2-1,n11+1:n22-1).*Ez(n1+1:n2-1,n11+1:n22-1)+(Hy(n1+1:n2-1,n11+1:n22-1)-Hy(n1:n2-2,n11+1:n22-1)-Hx(n1+1:n2-1,n11+1:n22-1)+Hx(n1+1:n2-1,n11:n22-2)).*D(n1+1:n2-1,n11+1:n22-1);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %forward boundary
    if n>=xdim-2-xsource
        Ez(xdim-2,3:1:ydim-3)=c0efffor*(Ez(xdim-3,3:1:ydim-3)+prev_prev_xfor(1,3:1:ydim-3))-prev_prev_x_minus_1for(1,3:1:ydim-3)+c2efffor*(prev_xfor(1,3:1:ydim-3)+prev_x_minus_1for(1,3:1:ydim-3))+c3efffor*(prev_x_minus_1for(1,2:1:ydim-4)+prev_x_minus_1for(1,4:1:ydim-2)+prev_xfor(1,2:1:ydim-4)+prev_xfor(1,4:1:ydim-2));
    end
    
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at forward boundary
    prev_prev_xfor=prev_xfor;
    prev_prev_x_minus_1for=prev_x_minus_1for;
    prev_xfor(1,1:1:ydim)=Ez(xdim-2,1:1:ydim);
    prev_x_minus_1for(1,1:1:ydim)=Ez(xdim-3,1:1:ydim);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %backward boundary
    if n>=xsource-3
        Ez(2,3:1:ydim-3)=-prev_prev_xrev(1,3:1:ydim-3)+c1effrev*(Ez(3,3:1:ydim-3)+prev_prev_x_minus_1rev(1,3:1:ydim-3))+c2effrev*(prev_xrev(1,3:1:ydim-3)+prev_x_minus_1rev(1,3:1:ydim-3))+c3effrev*(prev_x_minus_1rev(1,2:1:ydim-4)+prev_x_minus_1rev(1,4:1:ydim-2)+prev_xrev(1,2:1:ydim-4)+prev_xrev(1,4:1:ydim-2));
    end
    
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at backward boundary
    prev_prev_xrev=prev_xrev;
    prev_prev_x_minus_1rev=prev_x_minus_1rev;
    prev_xrev(1,1:1:ydim)=Ez(3,1:1:ydim);
    prev_x_minus_1rev(1,1:1:ydim)=Ez(2,1:1:ydim);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %upward boundary
    if n>=ydim-2-ysource
        Ez(3:1:xdim-3,ydim-2)=c0efffor*(Ez(3:1:xdim-3,ydim-3)+prev_prev_yfor(3:1:xdim-3,1))-prev_prev_y_minus_1for(3:1:xdim-3,1)+c2efffor*(prev_yfor(3:1:xdim-3,1)+prev_y_minus_1for(3:1:xdim-3,1))+c3efffor*(prev_y_minus_1for(2:1:xdim-4,1)+prev_y_minus_1for(4:1:xdim-2,1)+prev_yfor(2:1:xdim-4,1)+prev_yfor(4:1:xdim-2,1));
    end
        
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at upward boundary
    prev_prev_yfor=prev_yfor;
    prev_prev_y_minus_1for=prev_y_minus_1for;
    prev_yfor(1:1:xdim,1)=Ez(1:1:xdim,ydim-2);
    prev_y_minus_1for(1:1:xdim,1)=Ez(1:1:xdim,ydim-3);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %downward boundary
    if n>=ysource-3
        Ez(3:1:xdim-3,2)=-prev_prev_yrev(3:1:xdim-3,1)+c1effrev*(Ez(3:1:xdim-3,3)+prev_prev_y_minus_1rev(3:1:xdim-3,1))+c2effrev*(prev_yrev(3:1:xdim-3,1)+prev_y_minus_1rev(3:1:xdim-3,1))+c3effrev*(prev_y_minus_1rev(2:1:xdim-4,1)+prev_y_minus_1rev(4:1:xdim-2,1)+prev_yrev(2:1:xdim-4,1)+prev_yrev(4:1:xdim-2,1));
    end
       
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at downward boundary
    prev_prev_yrev=prev_yrev;
    prev_prev_y_minus_1rev=prev_y_minus_1rev;
    prev_yrev(1:1:xdim,1)=Ez(1:1:xdim,3);
    prev_y_minus_1rev(1:1:xdim,1)=Ez(1:1:xdim,2);
        
    %Mirroring of corner values taking the fact that corners are reached by the fields from the previous corners
    %in two time steps as S=1/sqrt(2) viz. sqrt(2)*delta(distance between two corners) is reached in 2 time steps
    Ez(2,2)=prev_prev_xrev(3);
    Ez(2,ydim-2)=prev_prev_xrev(ydim-3);
    Ez(xdim-2,2)=prev_prev_x_minus_1for(3);
    Ez(xdim-2,ydim-2)=prev_prev_x_minus_1for(ydim-3);
    
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
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with Mur absorbing boundary and at time = ',num2str(round(n*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
% [1] "Computational Electrodynamics - The Finite Difference Time Domain
%      Method" - Allen Taflove, Susan.c.Hagness, Third Edition, Artech House
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
