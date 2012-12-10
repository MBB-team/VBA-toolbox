tic
[T,R,Gdens,Gint] = propagator(0:0.2:8,0:0.2:16,2000,3,3);
toc
% close all;
dummy=10*log10(Gint/max(max(Gint)));
dbGint=nan(size(dummy));
dbGint(dummy>-30)=dummy(dummy>-30);
figure
pcolor(T,R,dbGint);
shading interp;
colorbar;
xlabel('t [ms]')
ylabel('r [mm]');