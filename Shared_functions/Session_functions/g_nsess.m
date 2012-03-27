function gx = g_nsess(x,P,u,in)

%in.dim_output
gx = zeros(in.dim_output,1);
%pause

for i=1:in.nsess
%     
% disp('ind.x')
% in.sess(i).ind.x
% disp('ind.phi')
% in.sess(i).ind.phi
% disp('ind.u')
% in.sess(i).ind.u
% disp('ind.gx')
% in.sess(i).ind.gx
% 
% disp('x')
% x
% disp('P')
% P
% disp('u')
% u
% 
% if isempty(in.sess(i).ind.phi)
% P = [];
% else
% P =P(in.sess(i).ind.phi);
% end


    gi = feval(...
        in.sess(i).g_fname,...
        x(in.sess(i).ind.x),...
        P(in.sess(i).ind.phi),...
        u(in.sess(i).ind.u),...
        in.sess(i).inG);
    
    
    gx(in.sess(i).ind.gx,1) = gi;

   % gx
    
end

