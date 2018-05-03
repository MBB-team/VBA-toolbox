function fx = f_HGFinGame(x,P,u,in)
% HGF-learner in a dyadic game
% function fx = f_HGFinGame(x,P,u,in)
% HGF stands for "Hierarchical Gaussian Filter" (cf. Mathys et al. 2014).
% In a dyadic game, HGF tracks P(o=1), i.e. the probaiblity that the
% opponent choses the option o=1 [see f_VBvolatile0.m].
% IN:
%   - x: the previous posterior sufficient statistics:
%   x(1)= o
%   x(2)= E[log-odds of P(o=1)]
%   x(3)= log V[log-odds of P(o=1)]
%   x(4)= E[log-volatility]
%   x(5)= log V[log-volatility]
%   - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
%   the notation of [Mathys et al. 2010].
%   - u: u(1) = o
%   - in: options set in options.inF
% OUT:
%   - fx: the updated posterior sufficient statistics (having observed o).


if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return
end
[fx] = f_VBvolatile0(x,P,u,in); % HGF update rule