function [kernels] = VBA_VolterraKernels(posterior,out,nt)
% Estimation of the system's 1st-order Volterra kernels
% function [my,vy,mg,vg,mx,vx] = VBA_VolterraKernels(posterior,out)
% IN:
%   - posterior,out: the output of the VBA system inversion
%   - nt: the maximul lag for the estimation of the Volterra kernel
% OUT:
%   - kernels: a structure, containing the fields .y, .g, .x, which refer
%   to the observed system, the simulated system observables, and the
%   simulated system hidden states. Each of these fields then contains the
%   following subfields:
%   - m: the coresponding estimated 1st-order Volterra kernels
%   - v: the estimation variance of 1st-order Volterra kernels
%   - R2: the percentage of variance explained, for each dimension

nu = out.dim.u;
n = out.dim.n;
p = out.dim.p;
if isequal(out.dim.n_t,1) || nu < 1 || ~any(out.u(:)) || n <1 || isempty(out.options.f_fname)
    % not a dynamical system
    kernels = [];
    return
end
if nargin <3
    try
        nt = min([options.kernelSize,16]);
    catch
        nt = min([out.dim.n_t,16]);
    end
end
np = nt*nu +1;

% deal with micro-resolution input
if out.options.microU && ~isequal(out.options.decim,1)
    u = zeros(nu,out.dim.n_t);
    for t=1:out.dim.n_t
        u(:,t) = mode(out.u(:,(t-1)*out.options.decim+1:t*out.options.decim),2);
    end
else % do not change input
    u = out.u;
end

% configurate kernel estimation
[opt.inG.dgdp] = VBA_conv2glm(u,nt); % build convolution matrices
opt.priors.muPhi = zeros(np,1);
opt.priors.SigmaPhi = 1e1*eye(np);
opt.verbose = 0;
opt.DisplayWin = 0;
% model dimensions
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = np;

% 1- Volterra kernels of observed data
kernels.y.m = zeros(p,nt,nu);
kernels.y.v = zeros(p,nt,nu);
kernels.y.R2 = zeros(p,1);
if out.options.verbose
    fprintf(1,['Deriving 1st-order Volterra kernels... '])
    fprintf(1,'%6.2f %%',0)
end
isbin = zeros(p,1);
for k = 1:p
    y = out.y(k,:)';
    for s=1:length(out.options.sources)
        if ismember(k,out.options.sources(s).out)
            isbin(k) = out.options.sources(s).type;
        end
    end
    if isbin(k)
        g_fname = @g_convSig;
        opt.binomial = 1;
    else
        g_fname = @g_conv0;
        opt.binomial = 0;
    end
    [pk,ok] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,opt);
    kernels.y.R2(k) = ok.fit.R2;
    if out.options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*k/(2*p+n)))
    end
    if out.options.DisplayWin
        try
            set(out.options.display.hm(2),'string',[num2str(floor(100*k/(2*p+n))),'%']);
            drawnow
        end
    end
    for i=1:nu;
        ind = (i-1)*nt+1:i*nt;
        kernels.y.m(k,:,i) = pk.muPhi(ind)';
        kernels.y.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
    end
end

% 2- Volterra kernels of simulated system
opt.binomial = 0;
kernels.g.m = zeros(p,nt,nu);
kernels.g.v = zeros(p,nt,nu);
kernels.g.R2 = zeros(p,1);
for k = 1:p
    if isbin(k)
        g_fname = @g_convSig;
    else
        g_fname = @g_conv0;
    end
    y = out.suffStat.gx(k,:)';
    [pk,ok] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,opt);
    kernels.g.R2(k) = ok.fit.R2;
    if out.options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*(k+p)/(2*p+n)))
    end
    if out.options.DisplayWin
        try
            set(out.options.display.hm(2),'string',[num2str(floor(100*(k+p)/(2*p+n))),'%']);
            drawnow
        end
    end
    for i=1:nu;
        ind = (i-1)*nt+1:i*nt;
        kernels.g.m(k,:,i) = pk.muPhi(ind)';
        kernels.g.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
    end
end

% 3- Volterra kernels of hidden states
g_fname = @g_conv0;
kernels.x.m = zeros(p,nt,nu);
kernels.x.v = zeros(p,nt,nu);
kernels.x.R2 = zeros(p,1);
for k = 1:n
    y = posterior.muX(k,:)';
    [pk,ok] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,opt);
    kernels.x.R2(k) = ok.fit.R2;
    if out.options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*(k+2*p)/(2*p+n)))
    end
    if out.options.DisplayWin
        try
            set(out.options.display.hm(2),'string',[num2str(floor(100*(k+2*p)/(2*p+n))),'%']);
            drawnow
        end
    end
    for i=1:nu;
        ind = (i-1)*nt+1:i*nt;
        kernels.x.m(k,:,i) = pk.muPhi(ind)';
        kernels.x.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
    end
end

if out.options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end


