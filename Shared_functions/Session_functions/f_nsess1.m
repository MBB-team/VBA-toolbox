function [fx] = f_nsess1(x,P,u,in)
%
tic
fx = zeros(size(x));
nsess = in.nsess;
nx_ps = size(x,1)/in.nsess;
nu_ps = size(u,1)/in.nsess;

for i=1:nsess

    [fi] = feval(...
        in.f_fname,... % the function to be called for session i
        x(nx_ps*(i-1)+1:nx_ps*i),... % the indices of the hidden states concerned by the evolution function for session i
        P(in.ind_theta(i,:)),... % the indices of the evolution parameters for session i
        u(nu_ps*(i-1)+1:nu_ps*i),... % the indices of the data for session i
        []); 
    
    fx(nx_ps*(i-1)+1:nx_ps*i) = fi;
        
end

toc