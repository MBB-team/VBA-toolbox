% FHN evolution function:
%
%                   [fx] = f_CaBBI_FHN(Xt,Theta,Isyn,inF)
%
% this function contains the equations and parameters of the FHN neuron
% model, and the over-time integration scheme


function [fx] = f_CaBBI_FHN(Xt,Theta,Isyn,inF)


dt = inF.dt;
Tau_Ca0 = inF.Tau_Ca0;


% Hidden variables
V = Xt(1);                                                                 % voltage
W = Xt(2);                                                                 % recovery variable
Ca = Xt(3);                                                                % [Ca2+] kinetics

% evolution parameters
Tau_Ca  = Tau_Ca0* exp(Theta(1));                                          % decay time-constant of [Ca2+] kinetics
Ca_base = Theta(2);                                                        % basal value of [Ca2+]
scale   = inF.k_Ca0* exp(Theta(3));                                        % scale parameter of [Ca2+] kinetis

% parameters of the evolution mode
g_Ca  = 5;                                                                 % maximal conductance of the calcium channel
E_Ca  = 120;                                                               % reversal potential of the calcium channel 

% compute the current of HVA calcium channel
V_scaled = (30*V) - 40;
s_inf = 1./( 1 + exp(- (V_scaled+45)/5 ) ); 
I_Ca  = g_Ca .* s_inf.* (V_scaled-E_Ca);                                   % HVA calcium current


% build the right hand side (RHS) of the evolution's differential equations
xdot = [ ( V - (power(V,3)/3) - W  + Isyn ) 
         (0.08*(V+0.7-(0.8*W)))
         ( -scale*I_Ca - ( (Ca-Ca_base)/Tau_Ca) ) 
         ];
    

fx = Xt + dt.*xdot;
 
