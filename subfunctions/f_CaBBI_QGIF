% QGIF evolution function:
%
%                   [fx] = f_CaBBI_QGIF(Xt,Theta,Isyn,inF)
%
% this function contains the equations and parameters of the QGIF neuron
% model, and the over-time integration scheme


function [fx] = f_CaBBI_QGIF(Xt,Theta,Isyn,inF)

dt = inF.dt;
Tau_Ca0 = inF.Tau_Ca0;

% hidden variables
V = Xt(1);
Ca = Xt(2);

% evolution parameters
Ca_base = Theta(2);                                              % basal value of [Ca2+]
scale = inF.k_Ca0*exp(Theta(3));                                 % scale parameter of [Ca2+] kinetics
Tau_Ca = Tau_Ca0* exp(Theta(1));                                 % decay time-constant of [Ca2+] kinetics


% other parameters of the evolution model
Cm = 1;                                                          % membrane capacitance
g_L = 0.05;                                                      % leak conductance
g_Ca  = 5;                                                       % maximal conductance of HVA calcium channel
E_Ca  = 120;                                                     % reversal potential of HVA calcium channel 
E_L = -65;                                                       % reversal potential of leak current
Vth = -59.9;                                                     % voltage threshold
Ith = 0.16;                                                      % threshold current
deltaV = 3.48;                                                   % spike slope factor


% Compute the currents
s_inf = 1./( 1 + exp(- (V+45)/5 ) );
I_Ca  = g_Ca .* s_inf.* (V-E_Ca);                                 % HVA calcium current
I_L = g_L * (V-E_L);                                              % leak current
f_sp = ((g_L/(2*deltaV))*power(V-Vth,2) ) + I_L -Ith;
V_peak = 30;                                
I_rep = 120;                                 
b_rep = 1e8;
c_rep = 0.1;  
a_rep =  1e12;
sigma_peak = 1; 
G_rep = a_rep*(1/(sqrt(pi)*sigma_peak))*exp(- power(V - V_peak,2)/power(sigma_peak,2) ); 
f_rep =I_rep* 1./ (1+exp(-b_rep*((G_rep)-c_rep))); 


% build the right hand side (RHS) of the evolution's differential equations
xdot = [ ((-I_L+ f_sp + Isyn) ./ Cm - (f_rep/dt) )
         ( -scale*I_Ca - ( (Ca-Ca_base)/Tau_Ca) )  
         ];


fx = Xt + dt.*xdot;
