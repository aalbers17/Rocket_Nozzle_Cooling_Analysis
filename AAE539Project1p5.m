% part v): Perform a local analysis at the throat to determine the
% temperatures Twg and Twl/ Here that the part a steady-state temperature
% distribution and use the techniques you would employ for a regenerative
% cooling jacket analysis. You will need to iterate on the heat transfer
% coefficients as afunction of wall temperature using correlations for the
% liquid and gas side. what Twg value do you obtain? The melting
% temperature of copper is 1350K. Do we have a problem?  What is the
% minimum required flow rate to avoid melting at the throat? 
clear all; clc; close all;
Dt = 0.004318000000000;
rho_water = 1000;
load A.mat At As1total
% Throat values from CEA: 
Mt = 1; % mach #
Pt = 7.9207;
rhot = 0.775;
A = pi*(.17/2).^2; % m^2
At = .17*.0254; % m^2
mdot = linspace(35,1050,30)*1000*.000063090; % kg/s
kcopper = 398; % W/m/K , Assume constant (changes very little w.r.t. temp)
kc = polyfit([linspace(0,100,11)],[.55575 .57864 .59803 .61450 .62856 .64060 .65091 .65969 .66702 .67288 .67703],3);
thw = .125*.0254; % thickness, m
% initialize temperature:
g = 9.81; 
muo =  9.851100000000001e-05;
Cpo = 4.4727e3;
Ko = .94197; % W/m/K
Pr_o = (Cpo*muo/Ko);
gammao = 1.144300000000000;
Pc = 200*6894.76; % Pa
Tc = 3.154570000000000e+03;
cstar = 1641.7; % m/s
Tt = zeros(10000,10);
Tt(1,:) = 298; 
Tl = zeros(10000,length(mdot));
Tl(1,:) = 2.990137362127620e+02; % Assume it has reached steady state temp, 299 K for 35GPM
Twl = zeros(10000,10);
Twl(1,:) = 298; 
Twg = zeros(10000,10);
Twg(1,:) = 1000;
muL = [];
mu_S = [];
hg = [];
sigma = [];
qdotg = [];
kwater = [];
Pr_water = []; Nud = []; Re_water = []; hl_seider = []; qdotl = []; checkmin = []; minEl = []; Twgfinal = [];
v_water = mdot./rho_water./(pi*Dt.^2/4); % m/s % constant, rho of water constant
Trg = 3.149399787489976e+03; % Tr at the throat: 
% q = hg(Tr-Twg) = k(Twg-Twl) = hl(Twl-Tl)
tol = 10000;
check = ones(10000,10)*10^8;
% check(:,1) = tol+1; % some arbitrary value larger than tolerance;
i = 1;
j = 1;
for j = 1:10
while check(i,j) >= tol & i <=9999
% Guess Twg:
Twg(i,j) = 1000 + .1*(i-1); % from trial and error, learned that 2000-3000K is close, begin there:
% Compute hg:
sigma(i,j) = 1./((.5.*Twg(i,j)/Tc.*(1+(gammao-1)./2.*Mt^2)+.5).^.68.*(1+(gammao-1)./2*Mt.^2).^.12);
hg = (0.026./(Dt.^.2).*(muo.^.2*Cpo./Pr_o.^.6)*(Pc*g./cstar).^.8).*sigma; % W/m^2/K
% solve for qdot:
qdotg(i,j) = hg(i,j)*(Trg-Twg(i,j));
% Compute Twl:
Twl(i,j) = Twg(i,j)-qdotg(i,j)*thw/kcopper;
% Compute hl with Tl(i-1): assume Tl(1) is 299K, what we found our Tavg to
% be for 35GPM: 
% viscosity:
muL = 1.74214e-3;
if Twl(i,j) <= 373
mu_S = exp(-3.7188 + 578.919/(-137.546+Twl(i,j)))/1000;
else
mu_S = 0.283636/1000;
end
% Cp of water:
Cp = (-203.606 + 1523.29*Tl(1)/1000 + -3196.413*(Tl(1)/1000)^2 + 2474.455*(Tl(1)/1000)^3 + 3.855326/(Tl(1)/1000)^2)/.01801488; % J/kg/K
% k of water: 
kwater(i,j) = kc(1)*Tl(1)^3+kc(2)*Tl(1)^2 + kc(3)*Tl(1) + kc(4);
% Pr of water:
Pr_water(i,j) = (Cp.*muL./kwater(i,j));
% Re of water: 
Re_water(i,j) = rho_water.*v_water(j).*Dt./muL; % not constant, need function of viscosity vs temp. 
% Seider-Tate:
Nud(i,j) = 0.027.*Re_water(i,j).^(4/5).*Pr_water(i,j).^(1/3).*(muL./mu_S).^.14;
hl_seider(i,j) = Nud(i,j).*kwater(i,j)./Dt; 
% Compute qdotl: 
qdotl(i,j) = hl_seider(i,j)*(Twl(i,j)-Tl(1));
% check to see if qdots are equal: 
check(i,j) = abs(qdotg(i,j)-qdotl(i,j));
i = i+1;
end
i = 1;
[checkmin(j) minEl(j)] = min(check(1:end-1,j));
Twgfinal(j) = Twg(minEl(j),end);
end


%  continue: 
% figure
plot((1:1:9999)./10+1000,check(1:end-1,1),'r*') % error plot to check progress
title('error plot for qdot(gas) = qdot(liq)')
ylim([0 2e5])
xlabel('Twg')
ylabel('qin-qout')
legend('35GPM')
figure
plot((1:1:9999)./10+1000,check(1:end-1,j),'r*') % error plot to check progress
title('error plot for qdot(gas) = qdot(liq)')
ylim([0 2e8])
xlabel('Twg')
ylabel('qin-qout')
legend('1050GPM')


% mdot required for Twg <= 1350K:  
Twgfinal1 = 1000+minEl(1)/10; % Kelvin at flow rate of 35GPM
Twgfinal2 = 1000+minEl(end)/10; % Kelvin at flow rate 240 GPM, this was the closest I was able to get to 1350K. 

