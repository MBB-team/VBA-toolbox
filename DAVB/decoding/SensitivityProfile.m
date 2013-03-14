function [profile,thetas_lbl]=SensitivityProfile(posterior,out,c)

%== initializations
theta = posterior.muTheta;
phi = posterior.muPhi;
x0=posterior.muX0;

options=out.options;
dim=options.dim;
inF=options.inF;
u=out.u(:,1:dim.n_t);

%== get canonical regressors
xBF.dt=out.options.inF.deltat*out.options.decim;
xBF.name='hrf';
xBF = spm_get_bf(xBF);
hrf = xBF.bf;
clear xBF;
% TODO: automatic
r{1} = (c(1,:)==1) - (c(1,:)==.01);
r{2} = (c(2,:)==1) - (c(2,:)==.01);
r{3} = r{1}.*r{2};

% identify source types
srcs = [];
for i=1:length(options.sources)
    srcs(options.sources(i).out) = options.sources(i).type;
end

%== loop over parmater
theta_i=[inF.indA inF.indB{:} inF.indC inF.indD{:}];
for t=theta_i
    
    % simul trajectory
    delta_theta = zeros(length(theta),1);
    delta_theta(t) = .05; %1e-2;
    y_p= simul(theta+theta.*delta_theta,x0,phi,options,u) ;
    y_m= simul(theta-theta.*delta_theta,x0,phi,options,u) ;

        % loop over observations

        for obs_i = 1:size(y_p,1)
              %for each regressor
                for r_i=1:length(r)
                if(srcs(obs_i)==0)
                    temp_r = conv(r{r_i},hrf'); % todo: no convolution if not fmRI data
                    temp_r = temp_r(1:length(r{r_i}));
                else
                    temp_r=r{r_i};
                end
                temp_r = [ones(size(temp_r))' temp_r'];
                b_p = regress(y_p(obs_i,:)',temp_r);
                b_m = regress(y_m(obs_i,:)',temp_r);
            
                profile_temp(r_i,obs_i) = (b_p(2)-b_m(2))/(2e-2);
            
        end
        end
    profile{t}=profile_temp;
end


inF= out.options.inF;
% theta_i=[inF.indA inF.indB{:} inF.indC inF.indD{:}];
% dydp_stim=dydp_stim(theta_i,:,:);
% 
thetas_lbl={};
for i=1:length(inF.indA)
    thetas_lbl{end+1} = sprintf('A_%d',i);
end
for k=1:length(inF.indB)
    for i=1:length(inF.indB{k})
    thetas_lbl{end+1} = sprintf('B^%d_%d',k,i);
    end
end
for i=1:length(inF.indC)
    thetas_lbl{end+1} = sprintf('C_%d',i);
end
for k=1:length(inF.indD)
    for i=1:length(inF.indD{k})
    thetas_lbl{end+1} = sprintf('D^%d_%d',k,i);
    end
end



    
    
end



function y=simul(theta,x0,phi,options,u)
%- initialization
dim=options.dim;
n_t=length(u);
x = zeros(dim.n,n_t);
y = zeros(dim.p,n_t);
% dgdp = zeros(dim.n_theta,dim.p,dim.n_t);

x(:,1) = x0;
y(:,1) = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim,1);

%- loop over time
for t = 2:n_t
        
    % Evaluate evolution function at past hidden state
    if dim.n > 0 
        [x(:,t),~,~] = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim,t);
    end
    
    % Evaluate observation function at current hidden state
    [y(:,t),~,~] = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim,t);

end

end

function drdt=myDiff(theta,x0,phi,options,u,regressors)

dt = 1e-3;
ref = getLandmarks(theta,x0,phi,options,u) ;

for t=1:length(theta)
    theta2 = theta;
    theta2(t) = theta2(t) + dt;
    pert = getLandmarks(theta2,x0,phi,options,u) ;
    dydp(t,:) = sum(pert-ref,2)'/dt;
end
end
    