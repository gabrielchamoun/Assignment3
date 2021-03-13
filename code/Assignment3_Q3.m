clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 1000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
dt = 10e-15;        % Time step 10fs -> -(Width/100)/vT
nt = 200;           % Simulation steps
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

% Bottleneck Settings
% Addition of box colliders to the sim
toggleDiffusive = 1; % 0 - Specular, 1 - Diffusive
% Providing opposing corners of desired boxes
box = [ % x1,y1, x2, y2 (Per box)
    0.8, 0, 1.2, 0.4;
    0.8, 0.6, 1.2, 1;
    ].*1e-7;

% Electric Field Settings
voltage = 0.8;  % in volts
nx = 200;       % Number of divisions in X
ny = 100;       % ... in Y
boxL = 40;      % length of box along X
boxW = 40;      % width of box along Y
sigma = 0.001;  % Conductivity in the box
V = Assignment2_Q2(nx, ny, boxL, boxW, sigma) * voltage;
vMap = reshape(V, [ny nx]);    % Reshaping Vector to a matrix
[Ex,Ey] = gradient(-vMap);
Ex = Ex ./ (Length/nx);    Ey = Ey ./ (Width/ny);    % Scaling E-Field

% Potential Surface Map
vMap = reshape(V, [ny nx]);    % Reshaping Vector to a matrix
figure('name', 'Voltage - FD Solution'), surf(vMap'); % Plotting

% Electric Field
figure('name', 'Electric Field');
quiver(Ex,Ey,1.1), title('Electric Field');


% Initializing all electrons with position and velocity
% Electron object stores position, velocity and drift
eObj = struct('x', 0, 'y', 0, 'vx', 0, 'vy', 0, ...
                        'vm', 0, 'vdx', 0, 'vdy', 0);
eCol = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eObj(i).x = rand()*Length;
    eObj(i).y = rand()*Width;
    for b = 1:size(box,1)
        while eObj(i).x <= max(box(b,1),box(b,3)) && ...
                eObj(i).x >= min(box(b,1),box(b,3)) && ...
                eObj(i).y <= max(box(b,2),box(b,4)) && ...
                eObj(i).y >= min(box(b,2),box(b,4))
            eObj(i).x = rand()*Length;
            eObj(i).y = rand()*Width;
        end
    end
    eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
    eObj(i).vm = sqrt(eObj(i).vx^2 + eObj(i).vy^2);
end

% Simulation loop
figure('name', 'MC Simulation');
t = 0;          % Init time
counter = 2;    % Init Sim Counter
while t < tStop
    t = t + dt;	% Incrementing Time
    
    % binning x and y coordinates
    [xpartbin, e] = discretize([eObj(:).x], nx);
    [ypartbin, e] = discretize([eObj(:).y], ny);
    
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
        Fx = qe * Ex(ypartbin(i), xpartbin(i));
        Fy = qe * Ey(ypartbin(i), xpartbin(i));
        eObj(i).vx = eObj(i).vx + 1/2*Fx/mn*dt;
        eObj(i).vy = eObj(i).vy + 1/2*Fy/mn*dt;
        
        % Drift velocity calculation
        dx = eObj(i).x - prevX(i);
        eObj(i).vdx = dx/dt;
        
        % Magnitude of velocity
        eObj(i).vm = sqrt(eObj(i).vx.^2 + eObj(i).vy.^2);
        
        % Box Collider Conditions.
        for b = 1:size(box,1)
            % Plotting the Boxes on the plot
            if counter == 2 % Only plot once
                boxX = [box(b,1), box(b,1), box(b,3), box(b,3), box(b,1)];
                boxY = [box(b,2), box(b,4), box(b,4), box(b,2), box(b,2)];
                fill(boxX, boxY, 'k');
            end

            % Boundaries defined by coordinates
            boundL = min(box(b,1),box(b,3));    % Left Wall
            boundR = max(box(b,1),box(b,3));    % Right Wall
            boundT = max(box(b,2),box(b,4));    % Top Wall
            boundB = min(box(b,2),box(b,4));    % Bottom Wall

            % Conditions set by boundaires
            condL = prevX(i) <= boundL && eObj(i).x >= boundL;
            condR = prevX(i) >= boundR && eObj(i).x <= boundR;
            condB = prevY(i) <= boundB && eObj(i).y >= boundB;
            condT = prevY(i) >= boundT && eObj(i).y <= boundT;
            % Corner Cases
            condTC = condT && condR || condT && condL;
            condBC = condB && condR || condB && condL;

            % Checking Left and Right of the box First. 
            % Start by confirming if within Top and Bottom Boundaries
            if eObj(i).y <= boundT && eObj(i).y >= boundB
                if condL || condR || condTC || condBC
                    if toggleDiffusive == 1 % If diffusion is toggled
                        posX = eObj(i).vx >= 0;	% Storing sign
                        posY = eObj(i).vy >= 0;
                        % re-thermalizing
                        eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
                        eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
                        % checking sign
                        if (eObj(i).vx >= 0) == posX
                            eObj(i).vx = -eObj(i).vx;
                        end
                        if (eObj(i).vy >= 0) ~= posY
                            eObj(i).vy = -eObj(i).vy;
                        end
                    else
                        eObj(i).vx = -eObj(i).vx;   % Just inverting
                    end
                    if condL    % Shifting reflection point
                        eObj(i).x = boundL - (eObj(i).x - boundL);
                    else
                        eObj(i).x = boundR + (boundR - eObj(i).x);
                    end
                end
            end
            % Checking Top and Bottom of the box. 
            % Start by confirming if within Left and Right Boundaries
            if eObj(i).x <= boundR && eObj(i).x >= boundL
                if condT || condB
                    if toggleDiffusive == 1
                        posX = eObj(i).vx >= 0;
                        posY = eObj(i).vy >= 0;  % Storing sign
                        % re-thermalizing
                        eObj(i).vx = (sqrt(vT^2 / 2)*randn(1,1));
                        eObj(i).vy = (sqrt(vT^2 / 2)*randn(1,1));
                        % checking sign
                        if (eObj(i).vx >= 0) ~= posX
                            eObj(i).vx = -eObj(i).vx;
                        end
                        if (eObj(i).vy >= 0) == posY
                            eObj(i).vy = -eObj(i).vy;
                        end
                    else
                        eObj(i).vy = -eObj(i).vy;
                    end
                    if condB    % Shifting reflection point
                        eObj(i).y = boundB - (eObj(i).y - boundB);
                    else
                        eObj(i).y = boundT + (boundT - eObj(i).y);
                    end
                end
            end
        end
        
        
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
    pause(0.001);	% Delay for animation
    
%     % Calculating drift current density over time
%     Jnx(:,counter) = qe*eCount*mean([eObj(:).vdx]);
%     
%     Time(:,counter) = t;        % Recording time in a vector
%     counter = counter + 1;      % Incrementing Sim Counter
end
hold off

% figure('name', 'Extracted Properties')
% 
% % Plotting drift current density over time
% subplot(2,2,3), plot(Time, Jnx, 'k', 'LineWidth',1.75);
% xlabel('Time (t)'), ylabel('Drift current density (A/m^2)');
% title('Average Current Behaviour');