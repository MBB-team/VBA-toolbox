function y = logit(x, a)
    if any(x >= a) || any(x <= 0)
        error('Error: argument out of range for function logit()');
    end
    
    y = log(x./(a-x));
return;