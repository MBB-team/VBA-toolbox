function [yu,dyudtheta] = u_GaussianBumps(Theta,u,inF)

% input function for free-form deterministic DCM

spread = exp(Theta(1));
scale = exp(Theta(inF.indscale));
centres = exp(Theta(inF.indcentres));
n = length(centres);
yu = 0;
dyudtheta = zeros(size(Theta));
for i=1:n
   dt = centres(i)-u;
   bump = scale(i)*exp(-0.5*dt.^2./spread);
   dyudtheta(1) = ...
       dyudtheta(1) + scale(i)*dt^2./(2*spread*exp(dt^2/(2*spread)));
   dyudtheta(inF.indscale(i)) = bump;
   dyudtheta(inF.indcentres(i)) = ...
       -scale(i)*centres(i)*dt./(spread*exp(dt^2/(2*spread)));
   yu = yu + bump;
end
