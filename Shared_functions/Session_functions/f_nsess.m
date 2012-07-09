function [fx] = f_nsess(x,P,u,in)
fx = zeros(size(x));
for i=1:in.nsess
    [fi] = feval(...
        in.sess(i).f_fname,... % the function to be called for session i
        x(in.sess(i).ind.x),... % the indices of the hidden states concerned by the evolution function for session i
        P(in.sess(i).ind.theta),... % the indices of the evolution parameters for session i
        u(in.sess(i).ind.u),... % the indices of the data for session i
        in.sess(i).inF); 
    fx(in.sess(i).ind.x) = fi; 
end

