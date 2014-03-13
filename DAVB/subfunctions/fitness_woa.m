function f = fitness_woa(x,P,u,in)
% obtain fitness for war of attrition players

n = size(x,1);

V = P(1); % value of the resource
if in.quaCost
    k = P(2); % scaling of quadratic time cost
else
    k = 1;
end

% get mean outcome for each pair of population
O = zeros(n,n);
m = k*in.gridt.^(in.quaCost+1);
for i=1:n
    mA = m(i);
    for j=1:n
        mB = m(j);
        if mA>mB
            O(i,j) = V - mB;
        elseif mA==mB
            O(i,j) = 0.5*V - mB;
        else
            O(i,j) = -mA;
        end
    end
end
% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = O*x;



