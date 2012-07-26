function gx = g_nsess(x,P,u,in)
gx = zeros(in.dim.p,1);
for i=1:in.nsess
    gi = feval(...
        in.sess(i).g_fname,...
        x(in.sess(i).ind.x),...
        P(in.sess(i).ind.phi),...
        u(in.sess(i).ind.u),...
        in.sess(i).inG);
    gx(in.sess(i).ind.gx,1) = gi;
end

% same but computing indices when possible instead of loading them
% nsess = in.nsess;
% nx_ps = size(x,1)/in.nsess;
% nu_ps = size(u,1)/in.nsess;
% for i=1:nsess
%     [gi] = feval(...
%         in.g_fname,... % the function to be called for session i
%         x(nx_ps*(i-1)+1:nx_ps*i),... % the indices of the hidden states concerned by the evolution function for session i
%         P(in.ind_phi(i,:)),... % the indices of the evolution parameters for session i
%         u(nu_ps*(i-1)+1:nu_ps*i),... % the indices of the data for session i
%         in.sess(i).inG);
%     gx(dim.p_ps*(i-1)+1:dim.p_ps*i) = gi;
% end
