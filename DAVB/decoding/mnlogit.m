function y = mnlogit(x,a,b,n)
  

    y=x;

    for t=1:size(x,2)
       ex=[exp(a*x(:,t)-b); n ]; 
       ytemp = (ex ./ (1+sum(ex))); 
       y(:,t)=ytemp(1:end-1);
    end

    

end