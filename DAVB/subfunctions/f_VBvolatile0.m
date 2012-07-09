function [fx] = f_VBvolatile0(x,P,u,in)
% computes Laplace-VB update rules for hidden states sufficient statistics
% [fx] = f_VBvolatile0(x,P,u,in)
% This is the one-step Markovian update rule for the posterior sufficient
% statistics of a volatile environment, as derived in [Mathys et al. 2010].
% Note: the state space of the response model also contains the previous
% posterior on x2, for predictions of the next input u.
% IN:
%   - x: the previous posterior sufficient statistics. These include, using
%   the notation of [Mathys et al. 2010]:
%   mu1 = x(1);
%   mu2 = x(2);
%   sa2 = x(3);
%   mu3 = x(4);
%   sa3 = x(5);
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: the current input to the learner.
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed u),
%   as well as book keeping of the current belief (c.f. response model).

% transform and define states and parameters
x(3) = exp(x(3));
x(5) = exp(x(5));
ka = in.lev2*sgm(P(1),in.kaub);
om = P(2);
th = sgm(P(3),in.thub);
rf = in.rf;

fx = zeros(size(x));

% 1st level
fx(1) = u(1);  

% 2nd level
s1h = sgm(x(2),1)*(1-sgm(x(2),1));
pe1 = fx(1) - sgm(x(2),1);
s2h = x(3) + exp(ka*x(4)+om);
fx(3) = 1/(s2h^-1 + s1h);
fx(2) = x(2) +fx(3)*pe1;

% 3rd level
pi3h = 1/(x(5)+th);
w2 = exp(ka*x(4)+om)/(x(3)+exp(ka*x(4)+om));
r2 = (exp(ka*x(4)+om)-x(3))/(exp(ka*x(4)+om)+x(3));
pe2 = (fx(3)+(fx(2)-x(2))^2)/(exp(ka*x(4)+om)+x(3)) -1;
fx(5) = 1/(pi3h + .5*ka^2*w2*(w2+r2*pe2));
if fx(5) <= 0
    if rf <= 0
        fx(4:5) = NaN;
    else
        fx(5) = 1/rf*x(5);
        fx(4) = x(4) + .5*ka*w2*fx(5)*pe2;
    end
else
        fx(4) = x(4) + .5*ka*w2*fx(5)*pe2;
end

% retransform states
fx(3) = log(fx(3));
fx(5) = log(fx(5));


