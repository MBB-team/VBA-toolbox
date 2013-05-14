function [y] = binomial_sample(y)


for t=1:size(y,2)
    temp = cumsum(y(:,t)) ;
    r= rand();
    ri = r>[0; temp] & r < [temp ; 1] ;
    y(:,t) = ri(1:end-1) ;
end
    

end