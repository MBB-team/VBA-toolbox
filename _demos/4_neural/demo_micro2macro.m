% This demo simulates the shape of the vmPFC fMRI response to value
% ratings, i.e. a truncated V (Mael Lebreton, personal com.).
% The underlying idea is as folows.
% Let us assume that subjects have a utility function defined over some
% continuous dimension x of the item they have to rate. This may be, for
% example, how sweet is the candy. This quantity is typically uncertain,
% i.e. it can be repsented by a (Gaussian) probability density function
% with mean mu and variance Sigma. According to the 'ensemble coding'
% hypothesis, this distribution is encoded via the input to vmPFC neurons,
% i.e. each vmPFC neuron receives a given value x=X, such that if one
% plotted the historgrams of the X's, one would recover p(x). Now let us
% assume that the firing of vmPFC neurons encodes the utility U(x=X) of the
% item's valued dimension. Thus, the average response of the neuron's
% ensemble is the expectation E[U(x)] of U(x) over p(x).
% The key trick is to use the typical concave utility functions of
% behavioural economics. This induces a non-trivial susceptibility to the
% uncertainty Sigma over x, such that E[U(x)] is a convex function of mu.

clear all
close all

a = 1;
b = 2;
dx = 2e-2;
x = -0.5:dx:2;

sx = VBA_sigmoid(x);
Ux = 2*sx-1;
d2Udx2 = 2*sx.*(1-sx).*(1-2.*sx);

% xp = x(x>=0);
% xn = x(x<0);
% Uxp = log(xp+1);
% Uxn = -log(1-xn);
% Ux = [Uxn,Uxp];
% d2Up = -1./(xp+1).^2;
% d2Un = 1./(1-xn).^2;
% d2Udx2 = [d2Un,d2Up];
% Ux = log(x+1);%1-exp(-a.*x);
% d2Udx2 = -1./(x+1).^2;%-a.^2*exp(-a.*x);

hf = figure('color',[1 1 1]);
ha = subplot(3,2,1,'parent',hf);
plot(ha,x,Ux), title(ha,'utility function')
hold(ha,'on')
plot([min(x),max(x)],[0,0],'r--')
box(ha,'off')
set(ha,'ygrid','on','xlim',[min(x),max(x)],'ylim',[min(Ux),max(Ux)])
xlabel(ha,'x')
ylabel(ha,'U(x)')


ha = subplot(3,2,2,'parent',hf);
plot(ha,x,d2Udx2), title(ha,'utility curvature')
hold(ha,'on')
plot([min(x),max(x)],[0,0],'r--')
box(ha,'off')
set(ha,'ygrid','on','xlim',[min(x),max(x)],'ylim',[min(Ux),max(Ux)])
xlabel(ha,'x')
ylabel(ha,'U''''(x)')


mu = x;
v = 5;%x.*(1-x);
EU = Ux + 0.5*v.*d2Udx2;
ha = subplot(3,2,3,'parent',hf);
plot(ha,mu,EU), title(ha,'expected utility')
hold(ha,'on')
plot([min(x),max(x)],[0,0],'r--')
box(ha,'off')
set(ha,'ygrid','on','xlim',[min(x),max(x)],'ylim',[min(Ux),max(Ux)])
xlabel(ha,'E[x]')
ylabel(ha,'E[U(x)]')

ha = subplot(3,2,4,'parent',hf);
plot(ha,Ux,EU), title(ha,'expected utility')
hold(ha,'on')
plot([min(x),max(x)],[0,0],'r--')
box(ha,'off')
set(ha,'ygrid','on','xlim',[min(Ux),max(Ux)],'ylim',[min(Ux),max(Ux)])
xlabel(ha,'U(E[x])')
ylabel(ha,'E[U(x)]')



N = 1e2;
EU2 = zeros(N,length(mu));
for i=1:length(mu)
    for ii=1:N
        X = mu(i) + sqrt(v).*randn;
        EU2(ii,i) = 2*VBA_sigmoid(X)-1;
    end
end
EU2 = mean(EU2,1);
ha = subplot(3,2,5,'parent',hf);
plot(ha,mu,EU2), title(ha,'E[U(x)] : mcmc')

ha = subplot(3,2,6,'parent',hf);
plot(ha,Ux,EU2), title(ha,'E[U(x)] vs U(E[x])')


VBA_getSubplots ();

