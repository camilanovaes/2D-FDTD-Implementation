clear all
close all

% Number of cells in Z direction
numCells = 200;
nt = 1500;

% Initialize fields to zero
Ex(numCells) = 0;
Hy(numCells) = 0;

% Gaussian pulse parameters
t0 = 40;     %Center of pulse
spread = 12; % Width of pulse

eps = 8.854e-12;
mi = 4*pi*1.e-7;
c = 2.99e8;

dx = 1.e-3;
dt = (dx/c)*0.99;
  
% Main loop (Loop C)
for t = 1:nt
  
  ti = t*dt;
  
  % Add E field source now at the center cell.
  % Excitation is a Gaussian pulse.
  Ex(numCells/2) = exp(-0.5*((t0-t)/spread)^2); % Hard source
 
  % (Loop A)
  % Calculate Ex field. Note that first 
  % cell is skipped in the loop, since we 
  % need to access k-1
  for k = 2:numCells
    Ex(k) = Ex(k) - 0.5*(Hy(k)-Hy(k-1));
  end

  % (Loop B)
  % Calculate Hy field. Note that last 
  % cell is skipped in the loop, since we 
  % need to access k+1
  for k = 1:numCells-1
    Hy(k) = Hy(k) - 0.5*(Ex(k+1)-(Ex(k)));
  end

  % Plot 
  plot(1:numCells, Ex+1,'r','linewidth',2);
  hold on 
  plot(1:numCells, Hy-1,'b','linewidth',2);
  axis([1 200 -2 2]);
  grid on;
  hold off
  legend('Ex','Hy');
  pause(0.001);
  
end
