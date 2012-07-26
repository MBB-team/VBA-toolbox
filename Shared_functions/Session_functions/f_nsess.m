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

% same but computing indices when possible instead of loading them
% nsess = in.nsess;
% nx_ps = size(x,1)/in.nsess;
% nu_ps = size(u,1)/in.nsess;
% for i=1:nsess
%     [fi] = feval(...
%         in.f_fname,... % the function to be called for session i
%         x(nx_ps*(i-1)+1:nx_ps*i),... % the indices of the hidden states concerned by the evolution function for session i
%         P(in.ind_theta(i,:)),... % the indices of the evolution parameters for session i
%         u(nu_ps*(i-1)+1:nu_ps*i),... % the indices of the data for session i
%         in.sess(i).inF); 
%     fx(nx_ps*(i-1)+1:nx_ps*i) = fi;
% end