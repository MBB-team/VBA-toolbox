function [posterior,suffStat] = VBA_VarParam(y,posterior,suffStat,dim,u,options)
% VB update of the precision hyperparameters
% function [posterior,suffStat] = VBA_VBVarParam(y,posterior,suffStat,dim,u,options)
%
% This function computes the natural parameters of the Gamma variational
% posterior pdf of the variance parameters (measurement noise and
% stochastic innovations).

if isequal(options.g_fname,@VBA_odeLim)
    clear VBA_odeLim
end

if options.DisplayWin
    set(options.display.hm(1),'string','VB: variance hyperparameters... ');
    set(options.display.hm(2),'string','0%');
    drawnow
end


%---- Measurement noise precision ----%
[gx,dG_dX,dG_dPhi] = VBA_evalFun('g',posterior.muX(:,1),posterior.muPhi,u(:,1),options,dim,1);

gsi = find([options.sources(:).type]==0);
for si=1:length(gsi)
    s_out = options.sources(gsi(si)).out ;
    i_s_out = find(options.isYout(s_out,1)==0);
    s_out = s_out(i_s_out);
    if ~isempty(s_out)
     % first store variance over predicted data
        iQyt=options.priors.iQy{1,si}(i_s_out,i_s_out);
    	ny = length(find(diag(iQyt)~=0));
        dy = y(s_out,1) - gx(s_out);
        dy2 = dy'*iQyt*dy; 
        posterior.a_sigma(si) = options.priors.a_sigma(si) + 0.5*ny;
        posterior.b_sigma(si) = options.priors.b_sigma(si) + 0.5*dy2 ; 
        if dim.n > 0
            posterior.b_sigma(si) = posterior.b_sigma(si) + 0.5*trace(dG_dX(:,s_out)*iQyt*dG_dX(:,s_out)'*posterior.SigmaX.current{1});
        end
        if dim.n_phi > 0
            posterior.b_sigma(si) = posterior.b_sigma(si) + 0.5*trace(dG_dPhi(:,s_out)*iQyt*dG_dPhi(:,s_out)'*posterior.SigmaPhi);
        end
    end
  
end

%- State noise precision
if dim.n >0
    a0 = posterior.a_alpha;
    b0 = posterior.b_alpha;
    iQx = VBA_inv(options.priors.iQx{1},options.params2update.x{1},'replace');
    nx = length(options.params2update.x{1});
    [fx,dF_dX,dF_dTheta] = VBA_evalFun('f',posterior.muX0,posterior.muTheta,u(:,1),options,dim,1);
    dx = posterior.muX(:,1) - fx;
    dx2 = dx'*iQx*dx;
    posterior.a_alpha = options.priors.a_alpha + 0.5*nx;
    posterior.b_alpha = options.priors.b_alpha + 0.5*dx2 + 0.5*trace(dF_dX*iQx*dF_dX'*posterior.SigmaX0) + 0.5*trace(iQx*posterior.SigmaX.current{1});
    if dim.n_theta >0
        posterior.b_alpha = posterior.b_alpha + 0.5*trace(dF_dTheta*iQx*dF_dTheta'*posterior.SigmaTheta);
    end
end

gsi = find([options.sources(:).type]==0);
for t=2:dim.n_t
[gx,dG_dX,dG_dPhi] = VBA_evalFun('g',posterior.muX(:,t),posterior.muPhi,u(:,t),options,dim,t);

%- Measurement noise precision         

for si=1:length(gsi)
    s_out = options.sources(gsi(si)).out ;
    is_s_out = find(options.isYout(s_out,t)==0);
    s_out = s_out(is_s_out);
    if ~isempty(s_out)
    % first store variance over predicted data
    iQyt=options.priors.iQy{t,si};
    iQyt = iQyt(is_s_out,is_s_out);
    ny = length(find(diag(iQyt)~=0));
    dy = y(s_out,t) - gx(s_out);
    dy2 = dy'*iQyt*dy; 
    posterior.a_sigma(si) = posterior.a_sigma(si) + 0.5*ny;
    posterior.b_sigma(si) = posterior.b_sigma(si) + 0.5*dy2; 
    if dim.n > 0
        posterior.b_sigma(si) = posterior.b_sigma(si) + 0.5*trace(dG_dX(:,s_out)*iQyt*dG_dX(:,s_out)'*posterior.SigmaX.current{t});
    end
    if dim.n_phi > 0
        posterior.b_sigma(si) = posterior.b_sigma(si) + 0.5*trace(dG_dPhi(:,s_out)*iQyt*dG_dPhi(:,s_out)'*posterior.SigmaPhi);
    end
    end
end
    
    %- State noise precision
    if dim.n>0 && t<dim.n_t
        iQx = VBA_inv(options.priors.iQx{t},options.params2update.x{t},'replace');
        nx = length(options.params2update.x{t});
        [fx,dF_dX,dF_dTheta] = VBA_evalFun('f',posterior.muX(:,t-1),posterior.muTheta,u(:,t),options,dim,t);
        dx = posterior.muX(:,t) - fx;
        dx2 = dx'*iQx*dx;
        posterior.a_alpha = posterior.a_alpha + 0.5*nx;
        posterior.b_alpha = posterior.b_alpha + 0.5*dx2 + 0.5*trace(dF_dX*iQx*dF_dX'*posterior.SigmaX.current{t-1}) + 0.5*trace(iQx*posterior.SigmaX.current{t}) - trace(iQx*dF_dX'*posterior.SigmaX.inter{t-1});
        if dim.n_theta >0
            posterior.b_alpha = posterior.b_alpha + 0.5*trace(dF_dTheta*iQx*dF_dTheta'*posterior.SigmaTheta);
        end
    end
    
    % Display progress
    if options.DisplayWin && mod(t,dim.n_t./10) < 1
        set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
        drawnow
    end
    
end

% regularize VB update
if dim.n>0 && posterior.b_alpha <=0
    posterior.a_alpha = a0;
    posterior.b_alpha = b0;
    VBA_disp('Warning: cancelling VB update of variance hyperparameter.',options);
end

if options.DisplayWin
    set(options.display.hm(2),'string','OK.');
    drawnow
end





