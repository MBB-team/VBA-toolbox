function [agx] = augmentgx(gx,in)

agx = gx;

if isfield(in,'diff')&& in.diff
    agx = [agx,[0;diff(gx)]];
end
if isfield(in,'int')&& in.int
    agx = [agx,cumsum(gx)];
end
if isfield(in,'div')&& in.div
    agx = [agx,gx./in.u];
end
if isfield(in,'log')&& in.log
    agx = [agx,zeros(size(agx))];
    agx(gx>0,end) = log(gx(gx>0)+1);
end