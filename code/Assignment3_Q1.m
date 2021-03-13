clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')

% MC Settings
eCount = 30000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
dt = 10e-15;        % Time step 10fs -> -(Width/100)/vT
nt = 100;           % Simulation steps
tStop = nt * dt;	% Stop Time

% Region settings
Length = 200e-9;
Width = 100e-9;

% Electron Properties
kB = 1.38066e-23;   % J/K 
m0 = 9.11e-31;      % Rest Mass
mn = 0.26*m0;       % Effective Mass
qe = 1.602e-19;     % Charge of an Electron 

% Thermal Velocity
Temp = 300;     % K
vT = sqrt((2*kB*Temp)/mn);

% Scattering Settings
% probability of scattering
toggleScatter = 1;    % Toggle Scattering OFF(0) ON(1)
tmin = 0.2e-12;
pScatter = 1 - exp(-dt/tmin);
scatterTracker = zeros(nt, eCount);

% Electric Field Settings
constV = 0.5;       % Constant voltage applied
Ex = constV/Length;  % electric field along the x-axis
Ey = 0;
Fx = qe * Ex;
Fy = qe * Ey;

% Initializing all electrons with a position and velocity
eObj = struct('x', 0, 'y', 0, 'vx', 0, 'vy', 0, ...
                        'vm', 0, 'vdx', 0, 'vdy', 0);
eCol = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eObj(i).x = rand()*Length;
    eObj(i).y = rand()*Width;
    eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vm = sqrt(eObj(i).vx^2 + eObj(i).vy^2);
end

t = 0;          % Init time
counter = 2;    % Init Counter
while t < tStop
    t = t + dt;	% Incrementing Time
    
    prevX = [eObj(:).x];    % Previous X location
    prevY = [eObj(:).y];    % Previous Y location
    for i = 1 : eCount
        
        % Updating position
        eObj(i).x = prevX(i) + eObj(i).vx * dt;
        eObj(i).y = prevY(i) + eObj(i).vy * dt;
        
        % Scattering effect - Randomize direction/magnitude of velocity
        % Probability of scattering based on p.
        if pScatter > rand() && toggleScatter	% 'if true'
            eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
            eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
            scatterTracker(counter,i) = t;
        else
            scatterTracker(counter,i) = scatterTracker(counter-1,i);
        end
        
        % Updating the velocities with the effect of E-Field
        eObj(i).vx = eObj(i).vx + 1/2*Fx/mn*dt;
        eObj(i).vy = eObj(i).vy + 1/2*Fy/mn*dt;
        
        % Drift velocity calculation
        dx = eObj(i).x - prevX(i);
        eObj(i).vdx = dx/dt;
        
        % Magnitude of velocity
        eObj(i).vm = sqrt(eObj(i).vx.^2 + eObj(i).vy.^2);
        
        % Top and Bottom Boundary Conditions
        if eObj(i).y > Width % y = 200nm boundary 
            diff = eObj(i).y - Width;
            eObj(i).y = Width - diff;
            eObj(i).vy = -eObj(i).vy;
        end
        if eObj(i).y < 0   % y = 0nm boundary
            diff = -eObj(i).y;
            eObj(i).y = diff;
            eObj(i).vy = -eObj(i).vy;
        end
        
        % Plotting Electrons
        if (i <= ePlotted)
            % Plotting previous and current position
            subplot(1,1,1);
            p = plot( [prevX(i), eObj(i).x], ...
                [prevY(i), eObj(i).y], ...
                'LineWidth', 1.5);
            p.Color = eCol(i,:);        % Attributing unique colour
            axis([0,Length,0,Width]);	% Plot Axis' set
            title('Electron Trajectories');
            hold on
        end
        
        % Left hand and Right hand side Boundary Conditions
        if eObj(i).x > Length   % x = 100nm boundary
            eObj(i).x = eObj(i).x - Length;
        end
        if eObj(i).x < 0   % x = 0nm boundary
            eObj(i).x = eObj(i).x + Length;
        end
        
    end
    pause(0.01);	% Delay for animation
    
    % Calculating the Average Temperature over time
    allTemperatures = ( ([eObj(:).vm].^2) .* mn ) ./ (kB*2);
    Temp(:,counter) = mean(allTemperatures);
    
    % Calculating drift current density over time
    Jnx(:,counter) = qe*eCount*mean([eObj(:).vdx]);
    
    Time(:,counter) = t;        % Recording time in a vector
    counter = counter + 1;      % Incrementing Sim Counter
end
hold off

figure('name', 'Extracted Properties')
divX = 50;  divY = 50;  % plot settings

% Plotting the Average Temperature over time
subplot(2,2,1), plot(Time, Temp, 'k', 'LineWidth',1.75);
xlabel('Time (t)'), ylabel('Temperature (K)');
title('Average Temperature');

% Plotting Temperature Map
xv = linspace(min([eObj(:).x]), max([eObj(:).x]), divX);
yv = linspace(min([eObj(:).y]), max([eObj(:).y]), divY);
[X,Y] = meshgrid(xv, yv);
Z = griddata([eObj(:).x],[eObj(:).y],allTemperatures,X,Y); % Mapping Temps
subplot(2, 2, 2);
surf(X, Y, Z),title('Temperature Map');
axis([0,Length,0,Width]);  % Plot Axis' set

% Plotting drift current density over time
subplot(2,2,3), plot(Time, Jnx, 'k', 'LineWidth',1.75);
xlabel('Time (t)'), ylabel('Drift current density (A/m^2)');
title('Average Current Behaviour');

% Plotting an Electron Density map
ptsX = linspace(0, Length, divX+1);
ptsY = linspace(0, Width, divY+1);
N = histcounts2([eObj(:).y], [eObj(:).x], ptsY, ptsX);  % Binning
subplot(2, 2, 4);
surf(X, Y, N),title('Electron Density Map');