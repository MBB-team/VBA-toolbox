function [fx] = f_Mathys_binary2(x,P,u,in)
% Compute the VB-laplace update rules for hidden states sufficient
% statistics as described in
% "A bayesian foundation for individual learning under uncertainty"
% IN:
%   - x: the previous posterior sufficient statistics. These include, using
%   the notation of [Mathys et al. 2010]:
%   mu1 = x(1);
%   mu2 = x(2);
%   sa2 = x(3);
%   mu3 = x(4);
%   sa3 = x(5);
%   Note: x(6) and x(7) are dummy states that serve in emitting a response
%   at each trial (c.f. response model) and do not play a role in the VB
%   update rules.
%       fx(6) = x(2);
%       fx(7) = x(3);
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: the current input to the learner.
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%   as well as book keeping of the current belief (c.f. response model).


% General remarks on notation
% - pe : prediction error
% - s : sigma
% - ka : kappa
% - om : omega
% - Xh : X 'hat'

fx = zeros(size(x));
fx = x;

%______________________________________________________________
%       FIRST DIMENSION
%______________________________________________________________

if (u(1) == 0) %action/choice 1

    
%________ update model 1 using feedback    
    
% transform and define states and parameters
x(3) = exp(x(3)); % sa2 : sigma2 needs to be positive -> exponential mapping
x(5) = exp(x(5)); % sa3 : sigma3 needs to be positive -> exponential mapping

ka = in.lev3.*sgm(P(1),in.kaub); % ka : kappa sigmoidal mapping + upper bound
om = P(2); % om : omega : no mapping
th = sgm(P(3),in.thub); % th : theta sigmoidal mapping + upper bound

rf = in.rf;


% 1st level
fx(1) = u(2); % [Eq. 21] % outcome

% 2nd level
s1h = sgm(x(2),1)*(1-sgm(x(2),1)); % [Eq. 26]
pe1 = fx(1) - sgm(x(2),1); % [Eq. 25] (through 21,24)
s2h = x(3) + exp(ka*x(4)+om); % [Eq. 27]
fx(3) = 1/(s2h^-1 + s1h); % [Eq. 22]
fx(2) = x(2) +fx(3)*pe1; % [Eq. 23]

% 3rd level
pi3h = 1/(x(5)+th); % [Eq. 31]
w2 = exp(ka*x(4)+om)/(x(3)+exp(ka*x(4)+om)); % [Eq. 32]
r2 = (exp(ka*x(4)+om)-x(3))/(exp(ka*x(4)+om)+x(3)); % [Eq. 33]
pe2 = (fx(3)+(fx(2)-x(2))^2)/(exp(ka*x(4)+om)+x(3)) -1; % [Eq. 34]
fx(5) = 1/(pi3h + .5*ka^2*w2*(w2+r2*pe2)); % [Eq. 29]


fx = fixVol(fx,x,w2,pe2,ka,th,rf,1);

% retransform states
x(3)  = log(x(3));
fx(3) = log(fx(3));
fx(5) = log(fx(5));

% keep track of previous reprensatation for prediction purposes
fx(6) = x(2);
fx(7) = x(3);



%___________ update model 2 with no feedback

% {
x(10) = exp(x(10)); % sa2 : sigma2 needs to be positive -> exponential mapping
x(12) = exp(x(12)); % sa3 : sigma3 needs to be positive -> exponential mapping
fx(10) = x(10) + exp(ka*x(11)+om);
fx(12) = x(12) + th;
fx(10) = log(fx(10));
fx(12) = log(fx(12));
%}
end


% TO DO:
% - include gradients dfdth, dfdx (TAKE TRANSPOSE: f is a row vector, -> option.checkGrads=1
%   for debugging)
%______________________________________________________________
%       SECOND DIMENSION
%______________________________________________________________
if (u(1) == 1) % action/choice 2


%________ update model 1 using feedback    

% transform and define states and parameters
x(10) = exp(x(10)); % sa2 : sigma2 needs to be positive -> exponential mapping
x(12) = exp(x(12)); % sa3 : sigma3 needs to be positive -> exponential mapping

ka = in.lev3.*sgm(P(1),in.kaub); % ka : kappa sigmoidal mapping + upper bound
om = P(2); % om : omega : no mapping
th = sgm(P(3),in.thub); % th : theta sigmoidal mapping + upper bound

rf = in.rf;

% 1st level
fx(8) = u(2); % [Eq. 21]

% 2nd level
s1h = sgm(x(9),1)*(1-sgm(x(9),1)); % [Eq. 26]
pe1 = fx(8) - sgm(x(9),1); % [Eq. 25] (through 21,24)
s2h = x(10) + exp(ka*x(11)+om); % [Eq. 27]
fx(10) = 1/(s2h^-1 + s1h); % [Eq. 22]
fx(9) = x(9) +fx(10)*pe1; % [Eq. 23]

% 3rd level
pi3h = 1/(x(12)+th); % [Eq. 31]
w2 = exp(ka*x(11)+om)/(x(10)+exp(ka*x(11)+om)); % [Eq. 32]
r2 = (exp(ka*x(11)+om)-x(10))/(exp(ka*x(11)+om)+x(10)); % [Eq. 33]
pe2 = (fx(10)+(fx(9)-x(9))^2)/(exp(ka*x(11)+om)+x(10)) -1; % [Eq. 34]
fx(12) = 1/(pi3h + .5*ka^2*w2*(w2+r2*pe2)); % [Eq. 29]

fx = fixVol(fx,x,w2,pe2,ka,th,rf,2);

% retransform states
x(10)  = log(x(10));
fx(10) = log(fx(10));
fx(12) = log(fx(12));

% keep track of previous reprensatation for prediction purposes
fx(13) = x(9);
fx(14) = x(10);

%___________ update model 1 with no feedback

% {
x(3) = exp(x(3)); % sa2 : sigma2 needs to be positive -> exponential mapping
x(5) = exp(x(5)); % sa3 : sigma3 needs to be positive -> exponential mapping
fx(3) = x(3) + exp(ka*x(4)+om);
fx(5) = x(5) + th;
fx(3) = log(fx(3));
fx(5) = log(fx(5));
%}




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

