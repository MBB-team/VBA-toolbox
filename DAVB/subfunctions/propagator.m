function [T,R,Gdens,Gint] = propagator(t_range,r_range,G_0,sigma,nubar)
%
[T,R]=meshgrid(t_range,r_range);
Gint=T*0;
Gdens=T*0;
delta_t=t_range(2)-t_range(1);
delta_r=r_range(2)-r_range(1);

for n=1:length(t_range)
    for m=1:length(r_range)
        t=T(m,n);
        r=R(m,n);
        Gint(m,n)=quad2d(@mypropagator,r,r+delta_r,t,t+delta_t);
        Gdens(m,n)=mypropagator(eps,eps);
    end
end


function G=mypropagator(rdr,tdt)
% tdt=t+dt;
% rdr=r+dr;
G_0 = 2000;
sigma = 3;
nubar = 3;
G=2*pi*G_0/(4*pi*sigma^2)./tdt.*exp(-((rdr.^2+nubar^2*tdt.^2)./(2*sigma*nubar*tdt)));
% G=2*pi*rdr*G_0/(4*pi*sigma^2)./tdt.*exp(-((rdr.^2+nubar^2*tdt.^2)./(2*sigma*nubar*tdt)));