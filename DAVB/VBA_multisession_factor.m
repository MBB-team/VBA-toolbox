function [posterior, out] = VBA_multisession_factor(posterior,out)

%% skip if no multisession
if ~(isfield(out.options,'multisession') && isfield(out.options.multisession,'expanded') && out.options.multisession.expanded)
    posteriors = posterior;
    return;
end

%% extract session_wise posteriors
Ts = [0 cumsum(out.options.multisession.split)];

multisession = out.options.inF{2};
out.options.inF = out.options.inF{1:end-1};
out.options.inG = out.options.inG{1:end-1};

dim = out.options.dim;

for i=1:numel(out.options.multisession.split)
    
    idx_t = (Ts(i)+1):Ts(i+1);   

    if dim.n > 0
        % state
        idx_X0 = multisession.indices.X0(:,i);
        
        posteriors(i).muX0 = posterior.muX0(idx_X0);
        posteriors(i).SigmaX0 = posterior.SigmaX0(idx_X0,idx_X0);
        
        posteriors(i).a_alpha =  posterior.a_alpha ;
        posteriors(i).b_alpha =  posterior.b_alpha ;
        
        posteriors(i).iQx =  posterior.iQx(:,idx_t); %cellfun(@(Q) Q(idx_X0,idx_X0),posterior.iQx(:,idx_t),'UniformOutput',false);
        posteriors(i).muX = posterior.muX(idx_X0,idx_t);
        posteriors(i).SigmaX.current =  cellfun(@(Q) Q(idx_X0,idx_X0),posterior.SigmaX.current(idx_t,:),'UniformOutput',false);
        posteriors(i).SigmaX.inter =  cellfun(@(Q) Q(idx_X0,idx_X0),posterior.SigmaX.inter(:,idx_t),'UniformOutput',false);
        
        % evolution
        idx_theta = multisession.indices.theta(:,i);
        
        posteriors(i).muTheta = posterior.muTheta(idx_theta);
        posteriors(i).SigmaTheta = posterior.SigmaTheta(idx_theta,idx_theta);
    end
        
    % observation
    idx_phi = multisession.indices.phi(:,i);
    
    posteriors(i).muPhi = posterior.muPhi(idx_phi);
    posteriors(i).SigmaPhi = posterior.SigmaPhi(idx_phi,idx_phi);
    
    % hyper param
    if isfield(posterior,'a_sigma')
        posteriors(i).a_sigma =  posterior.a_sigma ;
        posteriors(i).b_sigma =  posterior.b_sigma ;
    end
    
    if isfield(posterior,'a_alpha')
        posteriors(i).a_alpha =  posterior.a_alpha ;
        posteriors(i).b_alpha =  posterior.b_alpha ;
    end
    
    posteriors(i).iQy = posterior.iQy{idx_t,:};
    
end

posterior.perSession = posteriors;




end