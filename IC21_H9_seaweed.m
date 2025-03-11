%Station H9
%Hales et al 1994 Deep-Sea Research
clear all

Station= "Hales 1996 - H9";

%% definition of the spatial domain with two different resolutions
z_max=40e-2;     %[m] bottom sediment depth, should be a multiple of z_res
z_res=0.2e-2;     %[m] depth resolution
ndepths = 1 + z_max/z_res;     %[no unit] number of depth layers
depths = linspace(0, z_max, ndepths); %[m] depth
z_res = ones(size(depths))*z_res; %[m] depth resolution

%% definition of the temporal domain
stoptime = 10;       %[a] total timespan of the problem
interval=1/600000;          %[a] time resolution (1/60000 is nine minutes, 1/8760 is one hour; 1/365.2 is a day)
t_length=stoptime/interval;      %[no unit] number of time layers

%% bottom-water environmental conditions
T=2.2;      %[C] temperature
S=34.9;   %[psu] salinity
P=5312.4;    %[bar] pressure
rho_sw = gsw_rho(S,T,P);    %[kg/m^3] in situ seawater density computed from GSW toolbox

%% bottom-water values of dissolved species
dO2w=(266.6)*1e-6*rho_sw; %[mol/m3] dissolved oxygen from GLODAP at station location, bottom waters
dtalkw=(2342)*1e-6*rho_sw; %[mol/m3] dissolved oxygen from GLODAP at station location, bottom waters
dtCO2w=(2186)*1e-6*rho_sw; %[mol/m3] DIC from GLODAP at sation location, bottom waters
dtNO3w=(20.0668)*1e-6*rho_sw; %[mol/m3] nitrate from GLODAP at sation location, bottom waters
dtSO4w=(29264.2*S/35)*1e-6*rho_sw; %[mol/m3] computer from salinity (Millero, 2013)
dtPO4w=(1.3561)*1e-6*rho_sw; %[mol/m3] nitrate from GLODAP at sation location, bottom waters
dtNH4w=(0)*1e-6*rho_sw; %[mol/m3] assumed
dtH2Sw=(0)*1e-6*rho_sw; %[mol/m3] assumed
dFew=(0.5)*1e-9*rho_sw; %[mol/m3] typical for deep sea oxic bottom waters (Abadie et al., 2019)
dMnw=(0.5)*1e-9*rho_sw; %[mol/m3] typical for deep sea oxic bottom waters (Morton et al., 2019)
dtSiw=(120)*1e-6*rho_sw;  %[mol/m3] dissolved inorganic silica
dCaw=0.02128./40.087.*(S./1.80655)*rho_sw;  %[mol/m3] Ca, computed from salinity using Riley CG(1967)

%% depth-dependent porosity
phiBeta = 33;   %porosity attenuation coefficient
phiInf = 0.74;   %porosity at infinite depth
phi0 = 0.91;    %porosity at interface
phi = (phi0 - phiInf)*exp(-phiBeta*depths) + phiInf;   %porosity profile 
phiS=1-phi;   %solid volume fraction
tort=(1-2*log(phi)).^0.5;   %tortuosity from Boudreau (1996, GCA)
tort2=tort.^2;   %tortuosity squared

%% Redfield ratios
RC=106;     %Redfield ratio for C
RN=16;     %Redfield ratio for N 
RP=1;    % Redfield ratio for P
M_CH2O=30.031; %[g per mol]
M_NH3=17.031; %[g per mol]
M_H3PO4=97.994; %[g per mol]
M_OM=M_CH2O+(RN/RC)*M_NH3+(RP/RC)*M_H3PO4; %[g of OM per mol of OC] Organic Matter molar mass
%% Ratios for seaweed
RCs=630;     % C ratio for seaweed
RNs=70;     % N ratio for seaweed
RPs=1;    % P ratio for seaweed
M_CH2Os=30.031; %[g per mol]
M_NH3s=17.031; %[g per mol]
M_H3PO4s=97.994; %[g per mol]
M_OMs=M_CH2Os+(RNs/RCs)*M_NH3s+(RPs/RCs)*M_H3PO4s; %[g of OM per mol of OC] Organic Matter molar mass

%% solid fluxes and solid initial conditions
Foc=0.1790; %[mol/m2/a] flux of total organic carbon to the bottom 
Froc=Foc*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
Fsoc=Foc*0.27; %[mol/m2/a] flux of slow-decay organic carbon to the bottom 
Ffoc=Foc*0.70; %[mol/m2/a] flux of fast-decay organic carbon to the bottom 

FMnO2=0.0005; %typical for deep sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
FFeOH3=0.0005; %typical for deep sea oxic bottom waters (Archer et al., 2002; Boudreau, 1996)
Fcalcite=0.2; %[mol/m2/a] flux of calcite to the seafloor 
Faragonite=0; %[mol/m2/a] flux of aragonite to the seafloor
Fclay=26/360.31; %[mol/m2/a] flux of clay to the bottom: 360.31 is the molar mass of montmorillonite, typical deep sea clay

Ftot=Foc*M_OM+FMnO2*86.9368+FFeOH3*106.867+Fcalcite*100.0869+Faragonite*100.0869+Fclay*360.31; %[g/m2/a] total sediment flux 
v0=(Ftot)/(2.65e6*phiS(1));                                             %[m/a] bulk burial velocity at sediment-water interface
vinf=v0*phiS(1)/phiS(1,ndepths);                                    %[m/a]bulk burial velocity at the infinite depth
for j=1:ndepths
    u_nat(1,j)=vinf*phi(1,ndepths)/phi(1,j);                               %[m/a] porewater burial velocity
    w_nat(1,j)=vinf*phiS(1,ndepths)/phiS(1,j);                         %[m/a] solid burial velocity
end

%% solid fluxes and solid initial conditions for seaweed
Foc_s=4; %[mol/m2/a] flux of total organic carbon to the bottom 
Froc_s=Foc_s*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
Fsoc_s=Foc_s*0.27; %[mol/m2/a] flux of slow-decay organic carbon to the bottom 
Ffoc_s=Foc_s*0.70; %[mol/m2/a] flux of fast-decay organic carbon to the bottom

Ftots=Foc_s*M_OMs+FMnO2*86.9368+FFeOH3*106.867+Fcalcite*100.0869+Faragonite*100.0869+Fclay*360.31; %[g/m2/a] total sediment flux 
v0s=(Ftots)/(2.65e6*phiS(1));                                             %[m/a] bulk burial velocity at sediment-water interface
vinfs=v0s*phiS(1)/phiS(1,ndepths);                                    %[m/a]bulk burial velocity at the infinite depth
for j=1:ndepths
    u_s(1,j)=vinfs*phi(1,ndepths)/phi(1,j);                               %[m/a] porewater burial velocity
    w_s(1,j)=vinfs*phiS(1,ndepths)/phiS(1,j);                         %[m/a] solid burial velocity
end

%% reaction rates
kfast_0s = 70 ;
kslow_0s = 0.0001;
%kslowox =
%kslowanox =
%kfastox = 
%kfastanox =

%kslow_0=0.0001; %Sargassum SOURCE
%kslow_0 = 4.38; % aerobic at 4 degrees (Pedersen et al., 2021)
%kfast_0 = 70; %Sargassum
%kfast_0 = 44.9; % aerobic at 4 degrees Laminaria (Pedersen et al., 2021)

%% Start and stop time organic carbon influx seaweed
a = 0; %introduction of Organic carbon flux in years
b = 3/365;%stoptime; % end organic carbon flux in years
t_start = a / interval;
t_end = b / interval; 

%% diffusive boundary layer
dbl=0.938e-3;            %[m] thickness at location taken from Sulpis et al 2018 PNAS
load('ic_H9.mat');

rerun = 2;
time_saved_resolution=1/1000; %[a]