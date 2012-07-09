function  [ fx] = f_test_errors( x_t,P,u_t,in )
% Model : unidimensional massive object linked to a linear spring (fixed to zero, equilibrium length =0, constant = k) 
% and subject to a known external force
% random fluctuations affect both position and speed.
% INPUT 
% - x_t : speed, position
% - P : spring constant/mass, dissipation constant/mass
% - u_t : perturbation strength
% - in : structure
%   .dt : delta time
% OUTPUT
% - fx : speed, position


% Model dynamics
%-- update speed
% dx/dt = dx/dt(t-1) + a(t)*dt
% m*a(t) = -k*(x-x0) -f*v + u_t 
%(spring recall + dissipation + external force)
%-- update position
% x(t) = x(t-1) + dx/dt(t-1)*dt;

k = exp(P(1));
f = exp(P(2));

%----------- TESTING ERRORS

try
    e = in.error;
catch
    e = 0;
end

if e == 1 % out of bounds : parameter
    P(10);
elseif e == 2 % out of bounds : input
    u_t(10);
elseif e == 3 % out of bounds : hidden state
    x_t(10);
elseif e == 4 % erroneous calculation
    [1,1]+[1;1];
elseif e == 7 % nan input
    u_t = nan*u_t;
end


fx = x_t + [(-k*x_t(2) -f*x_t(1) + u_t);(x_t(1))]*in.dt;


if e == 5 % output nan
    fx = fx*nan;
elseif e == 6 % wrong dimensions
    fx = [fx;fx];
end
