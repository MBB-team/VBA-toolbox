function score = multinomial_accuracy(predict,data,isYout)
% MULTINOMIAL_ACCURACY
% score = multinomial_accuracy(predict,data)
% compute bias score of a multinomial classificator
% see Arzheimer & Evans (2013) A New Multinomial Accuracy Measure for Polling Bias, Political Analysis

if exist('isYout')
    isOutIdx = find(any(isYout));
    predict(:,isOutIdx) = [];
    data(:,isOutIdx)    = [];
end

B1 = compute_score(predict,data) ;

n = size(data,1) ;
chance = ones(n,1)/n;

B0 = compute_score(chance,data) ;

% score = exp(-(B1)) ;
% score = .5 + exp(-(B1)) - exp(-(B0));
% score = exp(-(B1/B0))/exp(-1) ;
% score = 1-(B1/B0) ;
score = 1-.5*(B1/B0) ;
%   exp(-(B0)) 


% score = 1-((1-exp(-(B1)))/(1-exp(-(B0)))) ;
% score = 1 - 1/(B1-B0) ;

end

function B = compute_score(p,d)
    p = mean(p,2) ;
    v = mean(d,2)    ;

    A = log( (p./(1-p)) ./ (v./(1-v))) ;

    B = nansum(v.*abs(A));
end