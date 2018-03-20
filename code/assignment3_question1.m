%%EMEKA PETERS - 100953293
%ELEC 4700 ASSIGNMENT 3 - Question 1
% clear all
clearvars
clearvars -GLOBAL
close all
format shorte


set(0, 'DefaultFigureWindowStyle', 'docked')
global C
%global Vx Vy x y Fx Fy AtomSpacing
%global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
%global LJEpsilon LJSigma Phi0 AtomType
%global MinX MaxX MinY MaxY PhiTot KETot
%global nAtoms0 nAtoms1 T T0 T1 MarkerSize
%global doPlotImage PlotCount map im PlotSize ScaleV ScaleF

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; % metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

effMass = 0.26 * C.m_0;

vTherm = sqrt((C.kb * 300) / effMass);

stdv = vTherm/(sqrt(2)); %Standard deviation for x and y velocities

dt = 7.5 * 10 ^ -15; % time-step value for iteration

temparr = zeros(1, 1000);
tmpx = (1:1:1000);

wid = 200 * 10 ^ -9; % x-boundaries
len = 200 * 10 ^ -9; % y-boundaries

volt = 0.1;
ex = volt/wid;
ey = 0;

size = 20000; %Number of Electrons
selectd = randi([1, size], 1, 60); %Selecting 60 random particles for plotting

% Assigning random positions to the electrons within the -
% boundaries
xPos = rand(1, size) .* wid;
yPos = rand(1, size) .* len;


isinbx = true;
while isinbx == true
   inbx = ((xPos <= (1.15 * wid/2) & (xPos >= (0.85 * wid/2))) & ((yPos < (len/3)) | yPos >= (2 * len/3)));
   if (sum(inbx) > 0)
       xPos(inbx) = rand(1, sum(inbx)) .* wid;
       yPos(inbx) = rand(1, sum(inbx)) .* len;
   else 
       isinbx = false;
       %break;
   end 
       
end

%Assigning a Random Velocity to the particles, following a
%Maxwell-Boltsmann distribution
velx = randn(1, size) .* stdv;
vely = randn(1, size) .* stdv;
vrms = sqrt((velx .^ 2) + (vely .^ 2));
vrmsarr = zeros(1, 1000);

%Accelerations
accx = zeros(1, size);
accy = zeros(1, size);

accx(1, :) = ex * (C.q_0/effMass);

driftvel = zeros(1, size);
electmob = zeros(1, size);

%Calculating the probability of scattering 
pscat = 1 - (exp((-1 * dt) / (0.2 * 10 ^ -12)));
tempr = 300;

%is2 = zeros(1, size);
noitr = 1000; % number of iterations
curvec = zeros(1, noitr);
timevec = zeros(1, noitr);
%curr = 0;
area = wid * len;
cardens = 10 ^ 15;


for i = 1:noitr
       
    %%This section models the electron motion 
    xPos(xPos >= wid) = xPos(xPos >= wid) - wid;
    xPos(xPos <= 0) = xPos(xPos <= 0) + wid;
    
    ylg = (yPos >= len);
    ylg1 = (yPos <= 0);
    
    vely(ylg) = -vely(ylg);
    vely(ylg1) = -vely(ylg1);
    
    xPosPrev = xPos;
    yPosPrev = yPos;
    
    velx = velx + (accx .* dt);
    
    xPos = xPosPrev + (velx .* dt);
    yPos = yPosPrev + (vely .* dt);
    
    
    vrms = sqrt((velx .^ 2) + (vely .^ 2));
    %vrmsarr(1, i) = vrms;
    tempr = (sqrt(2)*(mean(vrms) ^ 2) * effMass) / C.kb;
    temparr(1, i) = tempr;
    
    electmob = mean(vrms); % assumed electon mobility
    
    driftvel = electmob * ex; % drift velocity
    
    curr = area * cardens * (sum(driftvel)/size) * C.q_0;
    
    curvec(i) = curr;
    timevec(i) = i;
    
    %Scattering
    is = pscat > rand(1,size);
    
    velx(is) = randn .* stdv;
    vely(is) = randn .* stdv;

    
    figure (1);
    plot(xPos(selectd), yPos(selectd), '.');
    xlabel("x-Position");
    ylabel("y-Position");
    title(["Average Temperature = " num2str(tempr)]);
    
    xlim([0 wid]);
    ylim([0 len]);
    %pause(0.1);
    hold on
    
end

%Plot of Current vs time
figure (2)
plot(timevec, curvec);
xlabel('time');
ylabel('Average Current');
title('Average Current vs Time');
hold on


%Creating a plot of temperature versus timestep
figure (3)
plot(tmpx, temparr);
xlabel('time');
ylabel('Average Temperature');
title('Temperature vs Time');
hold on


[xgr, ygr] = meshgrid(0:(wid/30):wid, 0:(len/30):len);
elecmat = zeros(30, 30);
tempmat = zeros(30, 30);
numelec = 0;
totvel = 0;

for ii = 1:30
    xmin = xgr(1, ii);
    xmax = xgr(1, ii+1);
    for jj = 1:30
        ymin = ygr(jj, 1);
        ymax = ygr(jj+1, 1);
        for kk = 1:size
            if((xPos(kk) > xmin) && (xPos(kk) < xmax) && ((yPos(kk) > ymin) && yPos(kk) < ymax))
                numelec = numelec + 1;
                elecmat(ii, jj) = elecmat(ii, jj) + 1;
                totvel = totvel + sqrt((velx(kk) .^ 2) + (vely(kk) .^ 2));
                if(numelec ~= 0)
                    tempmat(ii, jj) = ((sqrt(2)*(totvel/numelec) ^ 2) * effMass) / C.kb;
                end
            end
        end
        totvel = 0;
        numelec = 0;
    end
end

%%Question 3.3 and 3.4 - Creating the plots of the histogram from question
%2, and the electron densitry map, along with the tempreature map.
figure(4); histogram(vrms, 10); title('Histogram of Thermal Velocities');
figure(5); surf(flipud(elecmat)); title('Electron Density Map');
figure(6); surf(flipud(tempmat)); title('Temperature Mat');

