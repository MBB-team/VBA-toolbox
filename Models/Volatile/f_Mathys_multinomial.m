function [fx] = f_Mathys_multinomial(x,P,u,in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel learning of multiple independent probabilities of binary outcome
% Probabilities are further used for decision making
% Here, inputs do not update every probability at each trial.
% A single probability is updated using the input at each trial
% (note that the others are updated)
% ----------
% Na : number of alternatives
% u : index of the chosen action (the only one updated using the feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the VB-laplace update rules for hidden states sufficient
% statistics as described in
% "A bayesian foundation for individual learning under uncertainty"
% IN:
%   - x: the previous posterior sufficient statistics. 
%   Multiple probabilities are tracked. 7 hidden states are associated to each.
%   The total number of hidden states is 7*Na.
%   For the first probability, using the notation of [Mathys et al. 2010],
%   the hiddens states are:
%   mu1 = x(1);
%   mu2 = x(2);
%   sa2 = x(3);
%   mu3 = x(4);
%   sa3 = x(5);
%   Note: x(6) and x(7) are dummy states that serve in emitting a response
%   at each trial (c.f. response model) and do not play a role in the VB
%   update rules.  fx(6) = x(2),  fx(7) = x(3);
%   The hidden states of all probabilities are concatenated.
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: 
%       1. the current input to the learner.
%       2. the index of the time series the input is related to. (starting
%       at 1)
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%   as well as book keeping of the current belief (c.f. response model).
% ----------
% General remarks on notation
% - pe : prediction error
% - s : sigma
% - ka : kappa
% - om : omega
% - Xh : X 'hat'

fx = zeros(size(x));fx = x; % init and copy previous states
Na = size(x,1)/7; % number of alternatives
ind = [1:7]+7*(u(1)-1); % indices concerned by updates by feedback

%________________________________________________    
%________ update model  using feedback________
%________________________________________________    

% transform and define states and parameters
x(ind(3)) = exp(x(ind(3))); % sa2 : sigma2 needs to be positive -> exponential mapping
x(ind(5)) = exp(x(ind(5))); % sa3 : sigma3 needs to be positive -> exponential mapping

ka = in.lev3.*sgm(P(1),in.kaub); % ka : kappa sigmoidal mapping + upper bound
om = P(2); % om : omega : no mapping
th = sgm(P(3),in.thub); % th : theta sigmoidal mapping + upper bound

rf = in.rf;


% 1st level
fx(ind(1)) =  u(1+u(1)) ; % [Eq. 21] % outcome

% 2nd level
s1h = sgm(x(ind(2)),1)*(1-sgm(x(ind(2)),1)); % [Eq. 26]
pe1 = fx(ind(1)) - sgm(x(ind(2)),1); % [Eq. 25] (through 21,24)
s2h = x(ind(3)) + exp(ka*x(ind(4))+om); % [Eq. 27]
fx(ind(3)) = 1/(s2h^-1 + s1h); % [Eq. 22]
fx(ind(2)) = x(ind(2)) +fx(ind(3))*pe1; % [Eq. 23]

% 3rd level
pi3h = 1/(x(ind(5))+th); % [Eq. 31]
w2 = exp(ka*x(ind(4))+om)/(x(ind(3))+exp(ka*x(ind(4))+om)); % [Eq. 32]
r2 = (exp(ka*x(ind(4))+om)-x(ind(3)))/(exp(ka*x(ind(4))+om)+x(ind(3))); % [Eq. 33]
pe2 = (fx(ind(3))+(fx(ind(2))-x(ind(2)))^2)/(exp(ka*x(ind(4))+om)+x(ind(3))) -1; % [Eq. 34]
fx(ind(5)) = 1/(pi3h + .5*ka^2*w2*(w2+r2*pe2)); % [Eq. 29]


fx = fixVol(fx,x,w2,pe2,ka,th,rf,u(1));

% retransform states
x(ind(3))  = log(x(ind(3)));
fx(ind(3)) = log(fx(ind(3)));
fx(ind(5)) = log(fx(ind(5)));

% keep track of previous reprensatation for prediction purposes
fx(ind(6)) = x(ind(2));
fx(ind(7)) = x(ind(3));


%________________________________________________    
%________ update other models without feedback___
%________________________________________________    
I = setdiff([1:Na],u(1)); % indices of the remaining models

for i = I
ind = [1:7]+7*(i-1);
x(ind(3)) = exp(x(ind(3))); % sa2 : sigma2 needs to be positive -> exponential mapping
x(ind(5)) = exp(x(ind(5))); % sa3 : sigma3 needs to be positive -> exponential mapping
fx(ind(3)) = x(ind(3)) + exp(ka*x(ind(4))+om);
fx(ind(5)) = x(ind(5)) + th;
fx(ind(3)) = log(fx(ind(3)));
fx(ind(5)) = log(fx(ind(5)));
end


function fx = fixVol(fx,x,w2,pe2,ka,th,rf,action)

% Extra treatment for the variance at level 3.
% This variance can't get too small.
% When update lead to a negative variance this update is discarded.
% Variance is reduced (as it also is with the discarded update) through a
% rescale. This is the Rescale Factor : rf
ind = [4:5] + (action-1)*7;

if fx(ind(2)) <= 0
    if rf <= 0
        fx(ind) = NaN;
    else
        fx(ind(2)) = x(ind(2)) + th;
        fx(ind(1)) = x(ind(1));% + (1/rf)*.5*ka*w2*fx(ind(2))*pe2; % [Eq. 30]
    end
else
        fx(ind(1)) = x(ind(1)) + .5*ka*w2*fx(ind(2))*pe2; % [Eq. 30]
end

