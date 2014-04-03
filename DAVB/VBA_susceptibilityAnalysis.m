function results = VBA_susceptibilityAnalysis_par(posterior,out,tapsOff,constraints)



% select only connexion of interest
% -------------------------------------------------------------------------
inF= out.options.inF;
Aself = inF.A;
Aself(Aself==1) = inF.indA;
true_As = setdiff(inF.indA,diag(Aself));
idxTheta=[true_As inF.indB{:} inF.indD{:}]; %inF.indC

thetas_lbl={};
for i=1:length(true_As)
    thetas_lbl{end+1} = sprintf('A_%d',i);
end
for k=1:length(inF.indB)
    for i=1:length(inF.indB{k})
    thetas_lbl{end+1} = sprintf('B^%d_%d',k,i);
    end
end
for k=1:length(inF.indD)
    for i=1:length(inF.indD{k})
    thetas_lbl{end+1} = sprintf('D^%d_%d',k,i);
    end
end

%% dimensions of interests
resp_idx = [out.options.sources(2:end).out];
nResps = numel(resp_idx);
nu = out.dim.u;
nw = numel(idxTheta);

%% flags
k_u = ones(nu,nu*nw);
for iu=1:nu
    k_u(iu,(iu-1)*nw+1:iu*nw) = 0 ;
end

k_w = ones(nw,nu*nw) ;
for iw=1:nw
    k_w(iw,iw+(0:nw:(nu-1)*nw)) = 0 ;
end

%% sub variances
out_temp = out;

% full model
ve1 = explainedVar(posterior,out,resp_idx,ones(nu,1),[]) ;

% random
ve0 = explainedVar(posterior,out,resp_idx,zeros(nu,1),[]) ;


% inputs
u_perms = flipud(full(spm_perm_mtx(nu)));

if strcmp(tapsOff,'oneatatime')
    u_perms = u_perms(sum(u_perms,2) == nu-1,:) ;
elseif strcmp(tapsOff,'allbutone')
    u_perms = u_perms(sum(u_perms,2) == 1,:) ;
elseif strcmp(tapsOff,'factorial')
    u_perms = u_perms(2:end-1,:) ;
end

parfor iu = 1:size(u_perms,1)
    v_u(iu,:) = explainedVar(posterior,out,resp_idx,u_perms(iu,:),[]) ;
end

% interactions
v_inter = [];
k_inter = [];
for iu = 1:size(u_perms,1)  
    parfor iw = 1:nw
        v_inter_iu(iw,:) = explainedVar(posterior,out,resp_idx,u_perms(iu,:),idxTheta(iw)) ; 
        k_temp = (1-k_w(iw,:))  .* (1-sum((1-k_u).*repmat((u_perms(iu,:)==0)',1,nu*nw)));
        k_inter_iu(iw,:) =  [k_temp (1-u_perms(iu,:))] ;
        
    end
    v_inter = [v_inter; v_inter_iu];
    k_inter = [k_inter; k_inter_iu]; 
end

% % parameters 
parfor iw = 1:nw
   v_w(iw,:) = explainedVar(posterior,out,resp_idx,ones(nu,1),idxTheta(iw)) ; 
end

%% inverting

f=k_inter;
f = [f ; (1-k_w) zeros(nw,nu)] ;
f(end+(1:size(u_perms,1)),nw*nu+(1:nu)) = 1-u_perms;

% f(:,end+1) = 1;

for i=1:nResps
    v_tot = ([v_inter(:,i); v_w(:,i); v_u(:,i)] - ve1(i) )/(ve0(i)-ve1(i)) ; %
%     v_tot = sqrt(v_tot);
    
    if constraints
        x_bar = fminsearch( @(x) norm(f*(x(:).^2) - v_tot), sqrt(ones(nw*nu+nu,1)/numel(nw*nu+nu)));
        x=x_bar(:).^2;
    else
        x = f\v_tot;
    end
    
    results.ev(:,i) = v_tot;
    results.ev_bar(:,i) = f*x;

%     x = x(1:end-1);

    contributions_w{i} = reshape(x(1:nw*nu),nw,nu);
    contributions_u{i} = x(nw*nu+(1:nu))';
    contributions_normoutput{i} = contributions_w{i} ./ repmat(sum(contributions_w{i}),nw,1) ;
    contributions_normoutput2{i} = contributions_w{i} ./ repmat(contributions_u{i},nw,1) ;

    contributions_normparam{i} = contributions_w{i} ./repmat(sum(contributions_w{i},2),1,nu) ;

    
    for iw=1:nw
        contributions_normparam{i}(iw,:) = contributions_w{i}(iw,:) ./ sum(contributions_w{i}(iw,:));
    end


end

%% results
    
results.out = out; 
results.posterior = posterior;

results.ev1 = ve1;
results.ev0 = ve0;

results.contributions_w = contributions_w;
results.contributions_u = contributions_u;
results.contributions_normoutput = contributions_normoutput;
results.contributions_normoutput2 = contributions_normoutput2;
results.contributions_normparam = contributions_normparam;

results.theta.idx = idxTheta;
results.theta.lbl = thetas_lbl;
results.f=f;





end

function v=explainedVar(posterior,out,resp_idx,u_switch,w_switch)
% v=0;
% return;

out.options.verbose = 0;
out.options.DisplayWin = 0;
out.options.inF.fast = true;

out.u = out.u .* repmat(u_switch(:),1,size(out.u,2));  
posterior.muTheta(w_switch) = 0;

% predict data
[yp,~,~,~,er] = simulateNLSS(...
    out.options.dim.n_t,...
    out.options.f_fname,...
    out.options.g_fname,...
    posterior.muTheta,...
    posterior.muPhi,...
    out.u,...
    Inf,...
    posterior.a_sigma./posterior.b_sigma,...
    out.options,...
    posterior.muX0);

g = yp-er;
y = out.y;% out.suffStat.gx; %

v = zeros(1,numel(resp_idx));
for i=1:numel(resp_idx)
    iY =resp_idx(i);
    in_idx = out.options.isYout(iY,:) == 0;
    v(i) = (var(y(iY,in_idx)-g(iY,in_idx)) / var(y(iY,in_idx) ))  ; 
end

end