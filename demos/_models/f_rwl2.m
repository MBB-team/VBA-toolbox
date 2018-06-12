function [fx] = f_rwl2(x,P,u,in)
% modified RW learning rule, with gains/losses asymmetry
% function [fx] = f_rwl2(x,P,u,in)
% This evolution function implements two different models (cf. in.model
% switch). If in.model = 'utility', different weights are given to
% negative, neutral and positive feedbacks. If in.model = 'learning',
% feedbacks are associated with different learning rates.
% IN:
%   - x: value of the 'go' option
%   - P: evolution parameters (depend upon the chosen model)
%   - u: feedback and last choice
%   - in: provides type of model and parameter indices
% OUT:
%   - fx: modified value (if last choice was 'go')

if u(in.inda) < 0.5 % last choice was 'nogo'
    fx = x; % no value update
elseif u(in.inda) >= 0.5 % last choice was 'go'
    switch in.model
        case 'utility' % weigths on feedbacks
            r = P(in.indR(u(in.indu)+2))*u(in.indu);
            alpha = VBA_sigmoid(P(in.indAlpha));
        case 'learning' % different learning rates
            r = P(in.indR)*u(in.indu);
            alpha = VBA_sigmoid(P(in.indAlpha(u(in.indu)+2)));
        case 'both'
            r = P(in.indR(u(in.indu)+2))*u(in.indu);
            alpha = VBA_sigmoid(P(in.indAlpha(u(in.indu)+2)));
        case 'none'
            r = P(in.indR)*u(in.indu);
            alpha = VBA_sigmoid(P(in.indAlpha));
        otherwise
            error
    end
    fx = x + alpha.*(r-x); % Rescorla-Wagner learning rule
end