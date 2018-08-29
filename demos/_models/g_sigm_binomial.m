function [Sx,dsdx,dsdp] = g_sigm_binomial(x,Phi,u,in)
% evaluates the sigmoid function for binomial data analysis

try
    in.x;
catch
    in.x = 0;
end

if in.x % for learning effects (sigmoid parameters evolve over time)
    [Sx, ~, dsdx] = VBA_sigmoid(u, ...
        'slope', exp(x(1)), ...
        'center', x(2), ...
        'derivatives', {'slope','center'});
    
    dsdx(1,:) = dsdx(1,:) * exp(x(1));  
    dsdp = [];
   
else
    [Sx, ~, dsdp] = VBA_sigmoid(u,...
        'slope', exp(Phi(1)), ...
        'center',Phi(2), ...
        'derivatives', {'slope','center'});
   
    dsdp(1,:) = dsdp(1,:) * exp(Phi(1)); 
    dsdx = [];
end


