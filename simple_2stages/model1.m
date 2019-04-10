%% Simple implementation of a stage structured biomass model
% Semi-chemostat phosphorus P, algal resource R,
% zooplankton juveniles J and adults A, 
% both feeding with type III functional response

%%
function dYdt = model1(t, Y)

d = 0.05;       %flow rate (1/d)
Pmax = 0.15;    %maximum P concentration (mg P/L)
q = 0.01;       %algal P content (mg P/mg C)
r = 1.5;        %maximum algal growth rate (1/d)
H = 0.0004;     %half-saturation constant (mg P/L)
aJ = 5;         %juvenile attack rate coefficient (L^b/((mg C)^b * d))
aA = 5;         %adult attack rate coefficient (L^b/((mg C)^b * d))
b = 1.5;        %Hill exponent (-)
hJ = 0.5;       %juvenile handling time (d)
hA = 0.5;       %adult handling time (d)
c = 0.6;        %conversion efficiency (-)
n = 0.01;       %maintenance losses (1/d)
z = 0.1;        %juvenile-to-adult body mass ratio (-)
mJ = 0.05;      %juvenile mortality rate (1/d)
mA = 0.05;      %adult mortality rate (1/d)

P = Y(1);       %ambient phosphorus concentration (mg P/L)
R = Y(2);       %algal density (mg C/L)
J = Y(3);       %juvenile density (mg C/L)
A = Y(4);       %adult density (mg C/L)

dPdt = d*(Pmax - P) - (q*r*P*R)/(H+P);
dRdt = (r*P*R)/(H+P) - (aJ*R^b*J)/(1+aJ*hJ*R^b) - ... 
    (aA*R^b*J)/(1+aA*hA*R^b) - d*R;
vJ = (c*aJ*R^b)/(1+aJ*hJ*R^b) - n;
vA = (c*aA*R^b)/(1+aA*hA*R^b) - n;
rJ = (vJ-mJ)/(1-z^(1-mJ/vJ));
dJdt = vJ*J + max(vA,0)*A - max(rJ,0)*J - mJ*J;
dAdt = max(rJ,0)*J + min(vA,0)*A - mA*A;

dYdt = [dPdt; dRdt; dJdt; dAdt];

end