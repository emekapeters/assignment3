%%EMEKA PETERS - 100953293
%ELEC 4700 ASSIGNMENT 3 - Questions 2 and 3
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



%%Question 2
%%Emeka Peters - 100953293
%ELEC 4700 - Assignment 2

numx = 80;
numy = 80;
Gmat = sparse((numx * numy), (numx * numy));
Vmat = zeros(1, (numx * numy));
vo = 0.8;

cond1 = 1; %Conductivity outside the boxes
cond2 = 1e-2; %Conductivity inside the boxes

box1dim = [(numx * 2/5), (numx * 3/5),  numy, (numy * 4/5)]; %Dimensions of box1 (bottom box)
box2dim = [(numx * 2/5), (numx * 3/5), 0, (numy * 1/5)]; %Dimensions of box2 (top box)
%{
Fractions were chosen for the box dimensions to give integer values after
multiplying with overall box size
%}



%Creating G-Matrix
for i = 1:numx

    for j = 1:numy
        
        n = j + (i - 1) * numy;
        
        if i == 1
            Gmat(n, :) = 0;
            Gmat(n, n) = 1;
            Vmat(1, n) = vo;
        elseif i == numx
            Gmat(n, :) = 0;
            Gmat(n, n) = 1;
            
        elseif j == 1 && i > 1 && i < numx
            
            if i == box1dim(1)
                Gmat(n, n) = -3;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond1;
            
            elseif i == box1dim(2)
                Gmat(n, n) = -3;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond2;
                
            elseif (i > box1dim(1) && i < box1dim(2))
                Gmat(n, n) = -3;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond2;
            else
                Gmat(n, n) = -3;
                Gmat(n, n + 1) = cond1;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond1;
            end
            
        elseif j == numy && i > 1 && i < numx
            
            if i == box1dim(1)
                Gmat(n, n) = -3;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond1;
            
            elseif i == box1dim(2)
                Gmat(n, n) = -3;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond2;
                
            elseif (i > box1dim(1) && i < box1dim(2)) 
                Gmat(n, n) = -3;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond2;
            else 
                Gmat(n, n) = -3;
                Gmat(n, n - 1) = cond1;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond1;
            end
            
        else
            
            if i == box1dim(1) && ((j < box2dim(4)) || (j > box1dim(4)))
                Gmat(n, n) = -4;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond1;
            
            elseif i == box1dim(2) && ((j < box2dim(4)) || (j > box1dim(4)))
                Gmat(n, n) = -4;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond2;
                
            elseif (i > box1dim(1) && i < box1dim(2) && ((j < box2dim(4)) || (j > box1dim(4))))
                Gmat(n, n) = -4;
                Gmat(n, n + 1) = cond2;
                Gmat(n, n - 1) = cond2;
                Gmat(n, n + numy) = cond2;
                Gmat(n, n - numy) = cond2;
            else
                Gmat(n, n) = -4;
                Gmat(n, n + 1) = cond1;
                Gmat(n, n - 1) = cond1;
                Gmat(n, n + numy) = cond1;
                Gmat(n, n - numy) = cond1;
            end
           
        end
    end
                

end 

figure (1);
spy(Gmat);
title('G Matrix Spy');
sol1 = Gmat\Vmat';
%surf(sol1);

actmat = zeros(numx, numy);



for i = 1:numx

    for j = 1:numy
        n = j + (i - 1) * numy;
        actmat(i, j) = sol1(n);
    end
end

figure (2);
surf(actmat);
title('Voltage Map with Bottlenecks');


[Efieldx, Efieldy] = gradient(actmat);
Efieldx = -Efieldx';
Efieldy = -Efieldy';

figure (3);
quiver(flipud(Efieldy), flipud(Efieldx));
title('Electric Field  Quiver Plot')
axis tight




%%Question 3

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

%accx(1, :) = ex * (C.q_0/effMass);

rnxPos = zeros(1, size);
rnyPos = zeros(1, size);

ind = zeros(1, size);

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

xPosPrev = xPos;
yPosPrev = yPos;

rnxPos = (floor(xPos .* 10^8)) + 1;
rnyPos = (floor(yPos .* 10^8)) + 1;
    
    %Indexing x and y Electric fields
ind = sub2ind([80 80], rnxPos, rnyPos);

for i = 1:noitr
       
    %%This section models the electron motion 
    xPos(xPos >= wid) = xPos(xPos >= wid) - wid;
    xPos(xPos <= 0) = xPos(xPos <= 0) + wid;
    
    ylg = (yPos >= len);
    ylg1 = (yPos <= 0);
    
    vely(ylg) = -vely(ylg);
    vely(ylg1) = -vely(ylg1);
    
    %box1dim = [(numx * 2/5), (numx * 3/5),  numy, (numy * 4/5)]; %Dimensions of box1 (bottom box)
    %box2dim = [(numx * 2/5), (numx * 3/5), 0, (numy * 1/5)]; %Dimensions of box2 (top box)
    
    inbx = ((xPos < (wid * 3/5) & (xPos > (wid * 2/5))) & ((yPos < (len/5)) | yPos > (4 * len/5)));
    inbx1 = (xPosPrev < (wid * 3/5) & (xPosPrev > (wid * 2/5)) & inbx);
    vely(inbx1) = -vely(inbx1);
    xPos(inbx1) = xPosPrev(inbx1);
    yPos(inbx1) = yPosPrev(inbx1);
    
    inbx2 = (xPosPrev > (wid * 3/5) & (xPosPrev < (wid * 2/5)) & inbx);
    velx(inbx2) = -velx(inbx2);
    xPos(inbx2) = xPosPrev(inbx2);
    yPos(inbx2) = yPosPrev(inbx2);  
    
    xPosPrev = xPos;
    yPosPrev = yPos;
    
    rnxPos = (floor(xPos(1, :) .* 10^8));
    rnyPos = (floor(yPos(1, :) .* 10^8));
    
    %Indexing x and y Electric fields
    %ind = sub2ind([numx numy], [rnxPos], [rnyPos]);
    
    %Solve for accelerations
    accx = 10^8 * Efieldy(ind) * (C.q_0/effMass);
    accy = 10^8 * Efieldx(ind) * (C.q_0/effMass);
    
    
    velx = velx + (accx .* dt);
    vely = vely + (accy .* dt);
    
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

    
    figure (4);
    plot(xPos(selectd), yPos(selectd), '.');
    xlabel("x-Position");
    ylabel("y-Position");
    title(["Average Temperature = " num2str(tempr)]);
    
    %box1dim = [(numx * 2/5), (numx * 3/5),  numy, (numy * 4/5)]; %Dimensions of box1 (bottom box)
    %box2dim = [(numx * 2/5), (numx * 3/5), 0, (numy * 1/5)]; %Dimensions of box2 (top box)
    
    line([wid * 2/5 wid * 2/5], [len 4*len/5]);
    line([wid * 3/5 wid * 3/5], [len 4*len/5]);
    line([wid * 2/5 wid * 3/5], [len len]);
    line([wid * 2/5 wid * 3/5], [4*len/5 4*len/5]);
    
    line([wid * 2/5 wid * 2/5], [0 len/5]);
    line([wid * 3/5 wid * 3/5], [0 len/5]);
    line([wid * 2/5 wid * 2/5], [0 0]);
    line([wid * 2/5 wid * 3/5], [len/5 len/5]); 
    
    xlim([0 wid]);
    ylim([0 len]);
    %pause(0.1);
    hold on
    
end

%Plot of Current vs time
figure (5)
plot(timevec, curvec);
xlabel('time');
ylabel('Average Current');
title('Average Current vs Time');
hold on


%Creating a plot of temperature versus timestep
figure (6)
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

%Electron Density and Temperature Maps
figure(7); histogram(vrms, 10); title('Histogram of Thermal Velocities');
figure(8); surf(flipud(elecmat)); title('Electron Density Map');
figure(9); surf(flipud(tempmat)); title('Temperature Mat');

