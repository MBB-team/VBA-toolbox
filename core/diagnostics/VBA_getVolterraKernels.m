function [kernels] = VBA_getVolterraKernels(posterior,out,nt)
% Estimation of the system's 1st-order Volterra kernels
% function [my,vy,mg,vg,mx,vx] = VBA_getVolterraKernels(posterior,out)
% IN:
%   - posterior,out: the output of the VBA system inversion. Note that the
%   Volterra kernels are estimated given the input stored in out.u.
%   - nt: the maximul lag for the estimation of the Volterra kernel
% OUT:
%   - kernels: a structure, containing the fields .y, .g, .x, which refer
%   to the observed system, the simulated system observables, and the
%   simulated system hidden states. Each of these fields then contains the
%   following subfields:
%   - m: the coresponding estimated 1st-order Volterra kernels
%   - v: the estimation variance of 1st-order Volterra kernels
%   - R2: the percentage of variance explained, for each dimension

nu = size(out.u,1);
n = out.dim.n;
p = out.dim.p;
if isequal(out.dim.n_t,1) || nu < 1 || ~any(out.u(:))
    % no time series or no input
    kernels = [];
    return
end
if nargin <3
    try
        nt = out.options.kernelSize;
    catch
        nt = min([out.dim.n_t,16]);
    end
end
np = nt*nu +1;

% deal with micro-resolution input
if out.options.microU && ~isequal(out.options.decim,1)
    u = zeros(nu,out.dim.n_t);
    for t=1:out.dim.n_t
        u(:,t) = mean(out.u(:,(t-1)*out.options.decim+1:t*out.options.decim),2);
    end
else % do not change input
    u = out.u(:,1:out.dim.n_t);
end
if VBA_isWeird (u)
    VBA_disp('Warning: zero-padding weird inputs for Volterra decompositions.',out.options)
    i0 = isinf(u) | isnan(u) | ~isreal(u);
    u(i0) = 0;
end
if isfield(out.options,'orthU') && out.options.orthU
    VBA_disp('Warning: orthogonalizing inputs for Volterra decompositions.',out.options)
    u = VBA_orth(u',0)';
end

% configurate kernel estimation
[opt.inG.dgdp] = VBA_conv2glm(u,nt); % build convolution matrices
if isfield(out.options,'detrendU') && ~~out.options.detrendU
    VBA_disp('Warning: detrending inputs for Volterra decompositions.',out.options)
    Trend = [];
    for i=0:out.options.detrendU
        Trend = [Trend,VBA_vec(1:out.dim.n_t).^i];
    end
    Trend = VBA_orth(Trend,1);
    opt.inG.dgdp = [opt.inG.dgdp;Trend(:,2:out.options.detrendU+1)'];
    np = np+out.options.detrendU;
end
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
for k = 1:p
    y = out.y(k,:)';
    isYout = out.options.isYout(k,:)';
    vy = var(y(isYout==0));
    if vy > eps % only if var(y)>0
        
        % find source tpye
        sInd = cellfun(@(x) ismember(k,x), {out.options.sources.out});
        switch out.options.sources(sInd).type
            case 0
                g_fname = @g_conv0;
                opt.sources.type = 0;
                opt.priors.a_sigma = 1;
                opt.priors.b_sigma = vy;
            case {1,2}
                g_fname = @g_convSig;
                opt.sources.type = 1;
        end
        
        opt.isYout = isYout;
        
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
        for i=1:nu
            ind = (i-1)*nt+1:i*nt;
            kernels.y.m(k,:,i) = pk.muPhi(ind)';
            kernels.y.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
        end
    end
end

% 2- Volterra kernels of simulated system
opt.sources.type = 0;
kernels.g.m = zeros(p,nt,nu);
kernels.g.v = zeros(p,nt,nu);
kernels.g.R2 = zeros(p,1);
for k = 1:p
    y = out.suffStat.gx(k,:)';
    if var(y)>eps % only if var(y)>0
        opt.priors.a_sigma = 1;
        opt.priors.b_sigma = var(y);
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
        for i=1:nu
            ind = (i-1)*nt+1:i*nt;
            kernels.g.m(k,:,i) = pk.muPhi(ind)';
            kernels.g.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
        end
    end
end

% 3- Volterra kernels of hidden states
if  n <1 || isempty(out.options.f_fname)
    kernels.x = [];
    if out.options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(' OK.')
        fprintf('\n')
    end
    return
end
g_fname = @g_conv0;
kernels.x.m = zeros(p,nt,nu);
kernels.x.v = zeros(p,nt,nu);
kernels.x.R2 = zeros(p,1);
for k = 1:n
    y = posterior.muX(k,:)';
    if var(y)>eps % only if var(y)>0
        opt.priors.a_sigma = 1;
        opt.priors.b_sigma = var(y);
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
        for i=1:nu
            ind = (i-1)*nt+1:i*nt;
            kernels.x.m(k,:,i) = pk.muPhi(ind)';
            kernels.x.v(k,:,i) = diag(pk.SigmaPhi(ind,ind))';
        end
    end
end

if out.options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end


