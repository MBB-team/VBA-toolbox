function [] = show_potential(posterior)
% evaluates double well potential
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
theta1 = posterior.muTheta(1);
theta2 = posterior.muTheta(2);


x = -4:1e-3:6;
y = (x-theta2).^2.*(x-theta1).^2;

hfig = findobj('tag','show');

if isempty(hfig)
    hfig = figure('tag','show');
else
%     clf(hfig)
end

figure(hfig)
ha = axes('parent',hfig);
plot(ha,x,y)
