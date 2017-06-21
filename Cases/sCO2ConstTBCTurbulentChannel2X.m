%%%%%Verification for Turbulent Channel Flow of sCO2,RANS/DES,const. ConstT BC %%
%%%%%%%%%%% vs. the Heat Transfer book by Nellis & Klein (pp. 635-668)%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Mahdi Nabil, June. 15th-16th, 2017 %%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
format compact

% Dimensions
W=750.02e-6;              % Channel Width  [m]
H=737.32e-6;              % Channel Height [m]
L=50e-3;                  % Channel Length [m]   
Lv=L*(2/5);               % Heated Portion of Channel Length used for Validation [m] 
AR=W/H;                   % Aspect Ratio (min to max dimension of the duct cross section, Alex confirmed)[-] (HT Book by Nellis & Klein, page 652)
Vol=W*H*L;                % Channel Volume [m^3]

% Hydraulic Diameter
A=W*H;                                  % Cross Sectional Area [m^2]
Perim=2*(W+H);                          % Wetted Perimeter [m]  
Dh=(4*A)/Perim;                         % Hydraulic Diameter [m]

% Flow Characteristics
Pr=4.875;                            % Prandtl Number [-]
mu=5.214e-5;                         % Dynamic Viscosity [kg/(m.s)] 
G=[731.196 733.02 516.016];          % Mass Flux [kg/(m^2.s)] (calc. in Paraview)  [ConstPropDES ConstPropmutDES ConstPropRANS]
rho=670;                             % Density [kg/(m^3)]  (const. in this simulation)
U=G./rho;                            % Flow Velocity [m/s] 
Re=(G.*Dh)./mu                       % Reynolds Number [-] (page 74 of Kyle's thesis)
L_plus=Lv./(Dh*Re);                  % Dimensionless Length [-] (HT Book by Nellis & Klein, page 651)
x_fd_h_T=Dh*(3.8*Re.^(1/6))*1000;    % Turbulent Hydrodynamic Entry Length [mm] (HT Book by Nellis & Klein, page 641)
x_fd_h_L=Dh*(0.06*Re)*1000;          % Laminar Hydrodynamic Entry Length [mm] (HT Book by Nellis & Klein, page 638)
x_fd_Th_T=Dh*7*1000;                 % Turbulent Thermal Entry Length [mm] (Sparrow & Cur, JHT, 1982, page 6/8)
x_fd_Th_L=Dh*(0.06*Re*Pr)*1000;      % Laminar Thermal Entry Length [mm] (HT Book by Nellis & Klein, page 647)

% Pressure Drop Calc.
f_fd_T=1./((0.79.*log(Re)-1.64).^2);                                                             % Friction factor for fully-developed turbulent flow (Smooth duct, e=0, eqn.5-63) [-]
f_fd_L=(96./Re).*(1-1.3553*AR+1.9467*(AR^2)-1.7012*(AR^3)+0.9564*(AR^4)-0.2537*(AR^5));          % Friction factor for fully-developed laminar flow (Smooth duct, e=0, eqn.5-60) [-]
f_avg_T=f_fd_T.*(1+(Dh/Lv).^0.7);                                                                % Apparent turbulent friction factor (eqn.5-66)[-]
f_avg_L=(4./Re).*((3.44./sqrt(L_plus))+(((1.25./(4.*L_plus))+((f_fd_L.*Re)./4)-(3.44./sqrt(L_plus)))./(1+(0.00021./(L_plus.^2)))));               % Apparent laminar friction factor (eqn.5-61)[-]
DP_T=(f_avg_T.*rho.*(U.^2).*Lv)/(2.*Dh)        % Theoretical Pressure drop along heated portion, Turbulent (HT Book by Nellis & Klein, page 660) [Pa] 
DP_L=(f_avg_L.*rho.*(U.^2).*Lv)/(2.*Dh)        % Theoretical Pressure drop along heated portion, Laminar (HT Book by Nellis & Klein, page 660) [Pa] 
DP_sim=[122.5 122.5 189.5]                           % Numerical Pressure drop along heated portion [Pa] (calc. in Paraview) [ConstPropDES ConstPropmutDES ConstPropRANS]
% Deviation from Correlation
AD_P=abs((DP_sim-DP_T)./DP_T)*100;             % Percentage of Absolute Deviation [-]

% Nu Correlations
Nu_G_corr_T=(((f_fd_T)./8).*(Re-1000).*Pr)./(1+12.7.*((Pr.^(2/3))-1).*sqrt((f_fd_T)./8))   % Local turbulent Nu [-] (eqn. 5-84 of Nellis, pp.667)
Nu_G_corr_avg=Nu_G_corr_T.*(1+(Lv/Dh).^-0.7);                                              % Average turbulent Nu [-] (eqn. 5-85 of Nellis, also Output of EES is the same)
Nu_Shah_corr_L=7.541*(1-2.610*AR+4.970*(AR^2)-5.119*(AR^3)+2.702*(AR^4)-0.548*(AR^5));    % Local laminar Nu, constT BC condition [-] (eqn. 5-82 of Nellis, pp.665)
Nu_DB_corr_T=0.023.*(Re.^(4/5)).*(Pr.^(0.4))             % Local turbulent Nu [-] (Dittus-Boelter equation)

% Nu Calc. from Simulation
Cp = 7221;             % Specific heat [J/(kg.K)], from  Fine2XConstProp case for T_f=313 K
k_f_avg = (Cp*mu)/Pr;     % Average thermal conductivity of fluid [W/(m.K)](which is almost constnt in our case since bulk fluid temp. is almost around 313 K)
q = [12929.508 12929.421 41600.161];             % Wall heat flux [W/m^2]  [ConstPropDES ConstPropmutDES ConstPropRANS]
T_f_avg = [304.042 304.042 305.99];      % Average fluid temp. [K]    [ConstPropDES ConstPropmutDES ConstPropRANS]
T_w_avg = [313.05 313.05 313.05];      % Average channel wall temp. [k] [ConstPropDES ConstPropmutDES ConstPropRANS]
Nu_sim  = (q.*Dh)./(k_f_avg.*(T_w_avg-T_f_avg))  % Nu from OF simulation [-]
% Deviation from Correlations
AD_Nu_G=abs((Nu_sim-Nu_G_corr_T)./Nu_G_corr_T)*100;           % Percentage of Absolute Deviation [-]
AD_Nu_DB=abs((Nu_sim-Nu_DB_corr_T)./Nu_DB_corr_T)*100;        % Percentage of Absolute Deviation [-]

