function  [ gx,dgdx,dgdP ] = g_ExpUtil( x_t,P,u_t,in )
% INPUT
% - x_t : two posterior moments of the mean and variance of u^(o)
% - P : inverse temperature
% - u_t : [useless]
% - in : []

beta = exp(P);
dQ = (x_t(1)-x_t(5));
gx = VBA_sigmoid( beta*dQ );

dgdx = zeros(8,1);
dgdx(1) = beta*gx*(1-gx);
dgdx(5) = -beta*gx*(1-gx);

dgdP = [beta*dQ*gx*(1-gx)];