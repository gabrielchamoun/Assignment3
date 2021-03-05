clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 1000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
dt = 10e-15;        % Time step 10fs -> -(Width/100)/vT
nt = 100;           % Simulation steps
tStop = nt * dt;	% Stop Time

% Region settings
Width = 200e-9;
Height = 100e-9;

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
mfp = vT * tmin;    % Mean Free Path
fprintf("tmin = %d m\n", tmin);
fprintf("Nominal Mean Free Path = %d m\n", mfp);

% Electric Field Settings
constV = 0.5;       % Constant voltage applied
Ex = constV/Width;  % electric field along the x-axis
Ey = 0;
Fx = qe * Ex;
Fy = qe * Ey;

% Initializing all electrons with a position and velocity
eObj = struct('x', 0, 'y', 0, 'vx', 0, 'vy', 0, 'vm', 0, 'vdx', 0);
eCol = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eObj(i).x = rand()*Width;
    eObj(i).y = rand()*Height;
    eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vm = sqrt(eObj(i).vx^2 + eObj(i).vy^2);
end

t = 0;          % Init time
counter = 2;    % Init Counter
while t < tStop
    xVect = 0;
    yVect = 0;
    t = t + dt;	% Incrementing Time
    for i = 1 : eCount
        
        % Updating position
        xVect(i) = eObj(i).x(counter-1) + eObj(i).vx * dt;
        yVect(i) = eObj(i).y(counter-1) + eObj(i).vy * dt;
        eObj(i).x(counter) = xVect(i);
        eObj(i).y(counter) = yVect(i);
        
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
        dx = eObj(i).x(counter) - eObj(i).x(counter-1);
        eObj(i).vdx = dx/dt;
        
        % Magnitude of velocity
        eObj(i).vm = sqrt(eObj(i).vx.^2 + eObj(i).vy.^2);
        
        % Top and Bottom Boundary Conditions
        if eObj(i).y(counter) > Height % y = 200nm boundary 
            diff = eObj(i).y(counter) - Height;
            eObj(i).y(counter) = Height - diff;
            eObj(i).vy = -eObj(i).vy;
        end
        if eObj(i).y(counter) < 0   % y = 0nm boundary
            diff = -eObj(i).y(counter);
            eObj(i).y(counter) = diff;
            eObj(i).vy = -eObj(i).vy;
        end
        
        % Plotting Electrons
        if (i <= ePlotted)
            % Plotting previous and current position
            subplot(3,1,1);
            p = plot( [eObj(i).x(counter-1), eObj(i).x(counter)], ...
                [eObj(i).y(counter-1), eObj(i).y(counter)], ...
                'LineWidth', 1.5);
            p.Color = eCol(i,:);        % Giving a unique colour to each
            axis([0,Width,0,Height]);	% Plot Axis' set
            title('Electron Trajectories');
            hold on
        end
        
        % Left hand and Right hand side Boundary Conditions
        if eObj(i).x(counter) > Width   % x = 100nm boundary
            eObj(i).x(counter) = eObj(i).x(counter) - Width;
        end
        if eObj(i).x(counter) < 0   % x = 0nm boundary
            eObj(i).x(counter) = eObj(i).x(counter) + Width;
        end
        
    end
    pause(0.01);	% Delay for animation
    
    
    % Plotting the Average Temperature over time
    Time(:,counter) = t;
    allTemperatures = ( ([eObj(:).vm].^2) .* mn ) ./ (kB*2);
    Temp(:,counter) = mean(allTemperatures); % Average Temperature
    subplot(3,2,3), plot(Time, Temp, 'k', 'LineWidth',1.75);
    xlabel('Time (t)'), ylabel('Temperature (K)');
    title('Average Temperature');
    
    % Plotting Temperature Map
    xv = linspace(min(xVect), max(xVect), 50);
    yv = linspace(min(yVect), max(yVect), 100);
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(xVect,yVect,allTemperatures,X,Y);  % Mapping Temps
    subplot(3, 2, 4);
    imagesc(xv,yv,Z),colorbar,title('Temperature Map');
    axis([0,Width,0,Height]);  % Plot Axis' set
    
    % Plotting drift current density over time
    Jnx(:,counter) = qe*eCount*mean([eObj(:).vdx]);
    subplot(3,2,5), plot(Time, Jnx, 'k', 'LineWidth',1.75);
    xlabel('Time (t)'), ylabel('Drift current density (A/m^2)');
    title('Average Current Behaviour');
    
    % Plotting an Electron Density map
    xPts = linspace(0, Width, 50);
    yPts = linspace(0, Height, 100);
    N = histcounts2(xVect, yVect, xPts, yPts);
    subplot(3, 2, 6);
    imagesc(xPts,yPts,N),colorbar,title('Electron Density Map');
    
    counter = counter + 1;      % Incrementing Sim Counter
end
hold off
