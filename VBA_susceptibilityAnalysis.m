function results = VBA_susceptibilityAnalysis(posterior,out,tapsOff)


% select only connexion of interest
% -------------------------------------------------------------------------
try 
    inF = out.options.inF{1};
catch
    inF = out.options.inF;
end

Aself = inF.A;
Aself(Aself==1) = inF.indA;
true_As = setdiff(inF.indA,diag(Aself));
idxTheta=[true_As inF.indB{:}  inF.indD{:}]; %inF.indC

thetas_lbl={};
for i=1:length(true_As)
    thetas_lbl{end+1} = sprintf('A_%d',i);
end
for k=1:length(inF.indB)
    for i=1:length(inF.indB{k})
    thetas_lbl{end+1} = sprintf('B^%d_%d',k,i);
    end
end
% for i=1:length(inF.indC)
%     thetas_lbl{end+1} = sprintf('C_%d',i);
% end
for k=1:length(inF.indD)
    for i=1:length(inF.indD{k})
    thetas_lbl{end+1} = sprintf('D^%d_%d',k,i);
    end
end

%% dimensions of interests
resp_idx = [out.options.sources(2:end).out];
nResps = numel(resp_idx);
u_idx = [1 2 5 6]; %1:size(out.u,1);% %  % % % %
nu = numel(u_idx);
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
ve1 = explainedVar(posterior,out,resp_idx,[],[]) ;

% random
ve0 = explainedVar(posterior,out,resp_idx,u_idx,[]) ;


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
    v_u(iu,:) = explainedVar(posterior,out,resp_idx,u_idx(u_perms(iu,:)==0),[]) ;
end

% interactions
v_inter = [];
k_inter = [];
for iu = 1:size(u_perms,1)  
    parfor iw = 1:nw
        v_inter_iu(iw,:) = explainedVar(posterior,out,resp_idx,u_idx(u_perms(iu,:)==0),idxTheta(iw)) ; 
        k_temp = (1-k_w(iw,:))  .* (1-sum((1-k_u).*repmat((u_perms(iu,:)==0)',1,nu*nw)));
        k_inter_iu(iw,:) =  [k_temp (1-u_perms(iu,:))] ;
        
    end
    v_inter = [v_inter; v_inter_iu];
    k_inter = [k_inter; k_inter_iu]; 
end

% % parameters 
parfor iw = 1:nw
   v_w(iw,:) = explainedVar(posterior,out,resp_idx,[],idxTheta(iw)) ; 
end

%% inverting

f=k_inter;
f = [f ; (1-k_w) zeros(nw,nu)] ;
f(end+(1:size(u_perms,1)),nw*nu+(1:nu)) = 1-u_perms;

% f(:,end+1) = 1;

for i=1:nResps
    v_tot = ([v_inter(:,i); v_w(:,i); v_u(:,i)] - ve1(i) )/(ve0(i)-ve1(i)) ; %

    x = f\v_tot;
        
    results.ev(:,i) = v_tot;
    results.ev_bar(:,i) = f*x;

%     x = x(1:end-1);

    contributions_w{i} = reshape(x(1:nw*nu),nw,nu);
    contributions_u{i} = x(nw*nu+(1:nu))';
    contributions_normoutput{i} = contributions_w{i} ./ repmat(sum(abs(contributions_w{i})),nw,1) ;
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

out.options.verbose = 0;
out.options.DisplayWin = 0;
inF.fast = true;

out.u(u_switch,:) = 0; 
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

g = yp-er; % predicted data under perturbation scheme
y = out.suffStat.gx; % original predicted data

v = zeros(1,numel(resp_idx));
for i=1:numel(resp_idx)
    iY =resp_idx(i);
    in_idx = find(out.options.isYout(iY,:) == 0);
%     v(i) = sqrt((var(y(iY,in_idx)-g(iY,in_idx)) / var(y(iY,in_idx) )))  ;
%       v(i) = mean(abs(y(iY,in_idx) - g(iY,in_idx))) /  mean(abs(y(iY,in_idx)));
      v(i) = mean(abs(y(iY,in_idx) - g(iY,in_idx))) ;
%       v(i) = std(y(iY,in_idx) - g(iY,in_idx)) / std(y(iY,in_idx));
end

end