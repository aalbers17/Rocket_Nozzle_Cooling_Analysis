% AAE539 Project 1: 
clear; close all; clc;

%% part i)
%Develop a code that calculates gas phase conditions (pressure,
% temperature, and density) as a function of axial distance along the
% nozzle. Create a nozzle contour that creates 100-200 points along the
% nozzle to provide significant resolution of heat flux gradients. 

% Procedure: run CEA for 5 ratios evenly distributed among both
% Sub/Supersonic regions. Then, polyfit for 200 points along axial distance
% to determine pressure, temperature, and density. 

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;
inp('type') = 'eq fr';              % Sets the type of CEA calculation
inp('p') = 200;                % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('o/f') = 3.5;               % Mixture ratio
inp('sup') = [1.412 2.824 4.236 5.6480 6.354 7.06];               % Supersonic area ratios
% inp('sub') = [5.292 4.704 3.528 2.352 1.176];
inp('fuel') = 'Paraffin_Wax';             % Fuel name from thermo.inp
inp('fuel_t') = 298;                % Fuel inlet temperature
inp('ox') = 'N2O4(L)';              % Ox name from thermo.inpj
inp('ox_t') = 298;                  % Ox inlet temperature
inp('file_name') = 'project1sup.inp';    % Input/output file name
if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').
data_eq = data('eq');
data_fr = data('fr');

% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.
temperature = squeeze(data_eq('t'));
pressure = squeeze(data_eq('p'));

% From CEA: 
xfit = [0 .075 .15 .3 .45 .6 .75 1 1.25 1.5 1.75 1.875 2.0];
Pfit = [13.789 13.688 13.661 13.559 13.255 11.083 7.9207 3.0997 1.0227 0.5669 0.37755 .32038 0.27688]*14.5038;
Tfit = [3154.57 3152.11 3151.45 3148.94 3141.39 3081.78 2969.85 2647.07 2236.09 2018.36 1875.41 1819.81 1771.51];
rhofit = [1.2571 1.249 1.2469 1.2387 1.2144 1.0387 0.775 0.3448 0.1356 0.0834 0.0598 .05229 0.0464];

% generate cubic spline to interpolate for 200 points: 
x = linspace(0,2,200);
dx = x(2)-x(1);
% generate radii from Solidworks sketch: 
dtheta = 3.29868987; 
theta = [];
rsubbend = [];
for i = 1:9
    theta(i) = 30.84893986-dtheta*(i);
    rsubbend(i) = .085+(.2-.2*cosd(theta(i)));
end
rsupbend = [];
theta = [];
for i = 1:7
    theta(i) = dtheta*(i);
    rsupbend(i) = .085+(.2-.2*cosd(theta(i)));
end
r = [linspace(.5,.1132954854,65) rsubbend .085 rsupbend linspace(.10099397,.6,118)];
D = r*2;
% r = [linspace(0.5,.17/2,75) linspace(.17/2+.0041,.6,124)]; 
P = spline(xfit,Pfit,x);
T = spline(xfit,Tfit,x);
rho = spline(xfit,rhofit,x);

% Tcst = polyfit(xfit,Tfit,3);
% rhocst = polyfit(xfit,rhofit,3);
% curve-fit with polynomial values and x: 
% P = pcst(1)*x.^3 + pcst(2)*x.^2 + pcst(3)*x + pcst(4) % *x.^7 + pcst(5)*x.^6 + pcst(6)*x.^5 + pcst(7)*x.^4+pcst(8)*x.^3 + pcst(9)*x.^2 + pcst(10)*x + pcst(11);
% T = Tcst(1)*x.^3 + Tcst(2)*x.^2 + Tcst(3)*x + Tcst(4) % *x.^7 + pcst(5)*x.^6 + pcst(6)*x.^5 + pcst(7)*x.^4+pcst(8)*x.^3 + pcst(9)*x.^2 + pcst(10)*x + pcst(11);
% rho = rhocst(1)*x.^3 + rhocst(2)*x.^2 + rhocst(3)*x + rhocst(4) % *x.^7 + pcst(5)*x.^6 + pcst(6)*x.^5 + pcst(7)*x.^4+pcst(8)*x.^3 + pcst(9)*x.^2 + pcst(10)*x + pcst(11);
% T = Tcst(1)*x.^10 + Tcst(2)*x.^9 + Tcst(3)*x.^8 + Tcst(4)*x.^7 + Tcst(5)*x.^6 + Tcst(6)*x.^5 + Tcst(7)*x.^4+Tcst(8)*x.^3 + Tcst(9)*x.^2 + Tcst(10)*x + Tcst(11);
% rho = rhocst(1)*x.^10 + rhocst(2)*x.^9 + rhocst(3)*x.^8 + rhocst(4)*x.^7 + rhocst(5)*x.^6 + rhocst(6)*x.^5 + rhocst(7)*x.^4+rhocst(8)*x.^3 + rhocst(9)*x.^2 + rhocst(10)*x + rhocst(11);

% legend('Pressure','Temperature','Density')


%% part ii) 
% modify your equation to solve the Bartz equation to predict the
% heat transfer coefficients, recovery temperature, and heat flux into the
% part for the assumed wall temperature of 1000K. 

% what we need: 
% Dt: throat diameter
% mu_o: chamber viscosity
% Cp_o: chamber specific heat
% Pr_o: chamber prandtl number
% Pc: chamber pressure
% cstar(x): C* at every location x
% g: gravity constant
% A(x): from geometry
% To: chamber temp
% gamma(x): gamma for every location x
% M(x): mach at every x
Twg = 1000; % assumed wall temp, K
g = 9.81; % m/s
Dt = .17*.0254; %m
muo = .000098511; % Pa*s
Cpo = 4472.7; % J/kg/K
Ko = .94197; % W/m/K
Pr_o = (Cpo*muo/Ko);
Pc = 200*6894.76; % Pa
Tc = 3154.57; % K
gammao = 1.1443; % from CEA
cstar = 1641.7; % m/s
% solve for A(x): 
A = pi*r.^2; % m^2
At = .17*.0254; % m^2
% solve for M(x) using CEA and spline interpolation: 
Mfit = [0 .113 .128 0.172 0.263 0.622 1 1.682 2.300 2.612 2.829 2.917 2.996];
afit = [1120.4  1119.8 1119.7 1119.2 1117.7 1105.7 1083.7 1025.7 953.6 910.3 879.3 866.7 855.5]; % speed of sound, needed to calculate v(x)

M = spline(xfit,Mfit,x); % Mach for every x:
a = spline(xfit,afit,x); % speed of sound for every x: 
% solve for v(x): 
v = M.*a; % m/s

% solve for hg(x): 
sigma = 1./((.5.*Twg/Tc.*(1+(gammao-1)./2.*M.^2)+.5).^.68.*(1+(gammao-1)./2*M.^2).^.12);
hg = (0.026./(Dt.^.2).*(muo.^.2*Cpo./Pr_o.^.6)*(Pc*g./cstar).^.8).*(At./A).^.9.*sigma; % W/m^2/K
hgmax = max(hg);

% solve for Tr: 
% need gamma(x), Pr(x), Cp(x), mu(x), K(x): same process as P,T,rho, M, a
Cpfit = [4.4727 4.4651   4.4630   4.4553   4.4318   4.2403 3.8580 2.7987   2.0553   1.9025   1.8495   1.8354   1.8253];
mufit = [0.000098511 0.000098458  0.000098444  0.000098391  0.000098229  0.000096950 0.000094522 0.000087298  0.000077588  0.000072218  0.000068600  0.000067172  0.000065920];
Kfit = [.94197 .93975   .93915   .93688   .93002   .87535  .77138 .49670   .28075   .22275   .19807   .19046   .18454];
gammafit = [1.1443 1.1443   1.1443   1.1444   1.1445   1.1457 1.1493 1.1702   1.2060   1.2188   1.2243   1.2259   1.2270];
Cp = spline(xfit,Cpfit,x);
mu = spline(xfit,mufit,x);
K = spline(xfit,Kfit,x);
gamma = spline(xfit,gammafit,x);
Pr_r = (Cp.*mu./K);
Re = rho.*v.*D./mu;
% because of turbulent free boundary layer: r = Pr^1/3
rf = Pr_r.^1/2;
Tr = T.*(1+(gamma-1)./2.*M.^2.*rf);
Trmax = max(Tr);
% solve for heat flux q: 
q = hg.*(Tr-Twg);
% plot results: 
figure
plot(x,r,x,r+.125)
title('Nozzle Profile')
xlabel('x [in]')
ylabel('r [in]')
legend('Gas-side','Water-side')
grid on
figure
plot(x,M,xfit,Mfit,'r*')
title('Mach vs. axial distance')
xlabel('x [in]')
ylabel('Mach')
legend('curve-fit','CEA points')
grid on
figure
plot(x,P,xfit,Pfit,'r*')
title('Pressure vs. axial distance')
xlabel('x [in]')
ylabel('P [psia]')
legend('curve-fit','CEA points')
grid on
figure
plot(x,T,xfit,Tfit,'r*')
title('Temperature vs. axial distance')
xlabel('x [in]')
ylabel('T [K]')
legend('curve-fit','CEA points')
grid on
figure 
plot(x,rho,xfit,rhofit,'r*')
title('Density vs. axial distance')
xlabel('x [in]')
ylabel('density [kg/m^3]')
legend('curve-fit','CEA points')
grid on
figure
plot(x,Tr)
title('Recovery Temp vs. axial distance')
xlabel('x [in]')
ylabel('Tr [K]')
% legend('curve-fit','CEA points')
grid on
figure
plot(x,hg)
title('h.t.c vs. axial distance')
xlabel('x [in]')
ylabel('hg [W/m^2/K]')
% legend('curve-fit','CEA points')
grid on
figure
plot(x,q)
title('Heat flux vs. axial distance')
xlabel('x [in]')
ylabel('q [W/m^2]')
% legend('curve-fit','CEA points')
grid on

%% part iii) 
% estimate an average heat transfer coefficient for the annulus
% perform an equilibrium balance with the heat going in vs the heat going
% out. then with everything known at the inlet of the annulus, solve for
% hl:
% step 1: obtain total heat out of system into annulus region - calculate
% and sum up all Qout for every point pertaining to annulus region. 
% step 2: solve for exit temperature of the water region using Q = mcpdT, where T1
% = 298K
% step 3: using Tavg (T1+T2/2) and Qtotal, solve for average hl

% Try again: Qnozzle = Qout 
% stored energy in nozzle: Qstored = mCpT
% mnozzle = 0.017133; % kg
% Cpnozzle = 390; % J/kg/K

% Qstored = mnozzle*Cpnozzle*Twg; % J

rho_water = 1000; % kg/m^3
% mdot = [5 10 35]*0.000063090*rho_water; % kg/s
% Cp = 4185.5; % J/kg/K for Tavg = 327, 312, 302 K. (little change w.r.t. temp)
Twl = 1000; % [K] assumed wall temp - liquid side
% annulus region can be found from x-points 26-150. find circumference for
% every point, multiply by distance of segment to next segment. 

for i = 26:149 % all but last segment to be dicretized (forward integration)
    As2(i-25) = pi*(D(i)+.25)*.0254*(sqrt(abs(x(i+1)*.0254-x(i)*.0254)^2+(abs(r(i)*.0254-r(i+1)*.0254))^2)); % m^2
%     As(i-25) = Asannulus(i-25) + pi*(1.25*.0254)*(x(i+1)*.0254-x(i)*.0254);
end
% Surface area of nozzle area: 
for i = 1:199 % all but last segment (forward step integration)
    As1(i) = pi*D(i)*.0254*(sqrt(abs(x(i+1)*.0254-x(i)*.0254)^2+(abs(r(i)*.0254-r(i+1)*.0254))^2)); % m^2
end
% As(1) = As(1) + (pi*((1.25*.0254)^2)/4-pi*((1.15*.0254)^2)/4); % for first discretization, add bottom region area
% As(end) = As(end) + (pi*((1.25*.0254)^2)/4 - pi*((1.18*.0254)^2)/4); % for last discretization, add top region area
As2total = sum(As2); % + pi*(1.25*.0254)*1.25*.0254 + (pi*(1.25^2)/4 - pi*(1.18^2)/4) + (pi*(1.25^2)/4-pi*(1.15^2)/4); % nozzle wall + top/bottom wall + outer wall area of water region
As1total = sum(As1);

Qin = q(1:end-1).*As1; % Qin to nozzle
Qtotal = sum(Qin); % total Watts
Qwaterin = linspace(Qtotal/200,Qtotal/200,100); % looking at half of system; As/2, Q/2...
% discretize water channel into 100 evenly sized chunks: 
Aswaterin = linspace(As2total/200, As2total/200,100); % m^2
Twater = zeros(100,100);
Cp = zeros(100,100);
Twater(1,:) = 298; % K 
% From NIST: (Shomate Equation)
Cp(1,:) = (-203.606 + 1523.29*Twater(1)/1000 + -3196.413*(Twater(1)/1000)^2 + 2474.455*(Twater(1)/1000)^3 + 3.855326/(Twater(1)/1000)^2)/.01801488; % J/kg/K
% assume mdot = 5 for now: 
mdot = linspace(0,35/2,100)*0.000063090*rho_water; % kg/s
for j = 1:100
    for i = 2:100
Twater(i,j) = Twater(i-1,j) + Qwaterin(1,i-1)/mdot(1,j)/Cp(i-1,j);
Cp(i,j) = (-203.606 + 1523.29*Twater(i,j)/1000 + -3196.413*(Twater(i,j)/1000)^2 + 2474.455*(Twater(i,j)/1000)^3 + 3.855326/(Twater(i,j)/1000)^2)/.01801488; % J/kg/K
    end
Tavg(j) = mean(Twater(:,j));
end

% % solve for hlavg:
hlavg = Qtotal/As2total./(Twl-Tavg); % W/m^2/K
% upon looking at these varying mdot values, we find that our hl values for
% htcs of reasonable  flow rates (Tfinal < 373K) are 14200-14000. 
% therefore, we will say that our hl value is: 14100. 
hl = 14100; 
% to validate this value, also calculate using Seider-Tate Correlation:
% define avg as location in which 50% of surface area is to the left and
% % right of the nozzle w.r.t. the axial direction. 
% mu_water = [5.04e-4 6.53e-4 8.145e-4]; % Pa-s, viscosity at Tavg
% mu_s = .038e-3; % Pa-s, viscosity of vapor on surface of nozzle at 1000K
% k_water = [.6406 .6286 .6145]; % W/m/K for Tavg of 327,312.302 K
% Pr_water = Cp.*mu_water/k_water; % prandtl number for various flow rates

% part iv): solve for bult temperature of the nozzle as a function of time:
% step 1: solve for Qstored: Qin-Qout = Qstored
% step 2: Solve for dT/dt: dT = (hgAs1(Tr-Twg) - hlAs2(Twl-Tl))/m/Cp
% step 3: dT(0) = 0, dT(10) = ? 
Tavg_final = [Tavg(15) Tavg(30) Tavg(100)];
mdot  = [5 10 35]*1000*.000063090;% grabbed from data corresponding to 5 10 35 GPM flows
% solve for Qin: 
Vnozzle = 4.269e-5; % m^3, from CAD
mnozzle = Vnozzle*8960; % V*rho 

% calculate nozzle with no flow, get volume of annulus: 
Vannulus = 3.665e-5; % m^3
mnozzlenoflow = mnozzle + Vannulus*8960; % kg
Cpnozzle = 390; % J/kg/K, assume constant

% Initialize integrated temperatures/heats: 
Twall = zeros(101,3);
Twallnoflow = zeros(101,1);
Twallnoflow = 298; % K
Twall(1,:) = 298; % initial temperature
Twater = zeros(101,3); % reinitialize Twater 
Twater(1,:) = 298; % Assume bulk water temp
Qout = zeros(101,3);
Qin = Qtotal; % no longer assuming equilibrium, however all heat generated is still this same value
% initialize Cp of water vs. temp again: 
Cp(1,:) = (-203.606 + 1523.29*Twater(1)/1000 + -3196.413*(Twater(1)/1000)^2 + 2474.455*(Twater(1)/1000)^3 + 3.855326/(Twater(1)/1000)^2)/.01801488; % J/kg/K
for i = 1:100 % from 0-10s: 100 points for better resolution
    for j = 1:3
Qout(i,j) = hl*As2total*(Twall(i,j)-Twater(i,j));
Twater(i+1,j) = Twater(i,j) + (1/10)*Qout(i,j)/(mdot(j)*(1/10))/Cp(i,j);
Twall(i+1,j) = Twall(i,j) + (1/10)*(Qin-Qout(i,j))/mnozzle/Cpnozzle;
Cp(i+1,j) = (-203.606 + 1523.29*Twater(i+1,j)/1000 + -3196.413*(Twater(i+1,j)/1000)^2 + 2474.455*(Twater(i+1,j)/1000)^3 + 3.855326/(Twater(i+1,j)/1000)^2)/.01801488; % J/kg/K
Twallnoflow(i+1,1) = Twallnoflow(i,1) + (1/10)*Qin/mnozzlenoflow/Cpnozzle;
    end
end
figure
plot(0:.1:10,Twall(:,1),0:.1:10,Twall(:,2),0:.1:10,Twall(:,3))
xlim([0 10]);
ylim([0 1200]);
title('Nozzle temperature vs. time')
xlabel('Time [s]')
ylabel('Temperature [K]')
legend('5 GPM','10 GPM','35 GPM')
grid on
figure 
plot(0:.1:10,Twallnoflow(:,1));
title('Nozzle with no cooling Chamber vs. time')
xlabel('Time [s]');
ylabel('Temperature [K]');
grid on

% part v)on separate script
 

