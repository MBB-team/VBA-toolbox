function [x,y,xhat,vx,yhat,vy] = VBA_comparePredictions(n_t,theta,phi,u,alpha,sigma,options,posterior,dim)
% compares inferred predictive density and real samples

if isempty(u)
    if options.microU
        u = zeros(1,n_t*options.decim);
    else
        u = zeros(1,n_t);
    end
end

% first get real system's sample path
f_fname = options.f_fname;
g_fname = options.g_fname;
opt = options;
opt.priors.muX0 = posterior.muX(:,end);
[y,x] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,opt,opt.priors.muX0);

% then get prediction from posterior onwards
post = posterior;
post.muX0 = posterior.muX(:,end);
post.SigmaX0 = posterior.SigmaX.current{end};
[xhat,vx,suffStat] = VBA_EKF(y,u,post,dim,opt,2);
yhat = suffStat.gx;
vy = suffStat.vy;

%display predictions
hf = figure(...
    'name','Predictive density: hidden states',...
    'color',[1 1 1],...
    'menubar','none');
hs = axes('parent',hf,'nextplot','add');
[haf,hf] = plotUncertainTimeSeries(xhat,getVar(vx,length(vx)),[],hs);
plot(x','--')
title('predicted hidden states')

hf = figure(...
    'name','Predictive density: observations',...
    'color',[1 1 1],...
    'menubar','none');
hs = axes('parent',hf,'nextplot','add');
[haf,hf] = plotUncertainTimeSeries(yhat,vy,[],hs);
plot(y','--')
title('predicted measurements')


function [] = dox(e1,e2)
ud = get(e1,'userdata');
scale = str2double(get(e1,'string'));
mu = ud.mu;
SX = ud.var;
n = size(mu,1);
for i = 1:n
    yp = [mu(i,:)+scale.*sqrt(SX(i,:)),fliplr(mu(i,:)-scale.*sqrt(SX(i,:)))];
    set(ud.hf(i),'ydata',yp)
end

function V = getVar(Sigma,indEnd)
if iscell(Sigma)
    n = size(Sigma{1},1);
    V = zeros(n,indEnd);
    for t=1:indEnd
        V(:,t) = diag(Sigma{t});
    end
else
    V = diag(Sigma);
end