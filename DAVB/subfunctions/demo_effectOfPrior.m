% Demo for Savage-Dickey nested model comparison

close all
clear variables

% Choose basic settings for simulations
sigma = 1e0;            % precision 
g_fname = @g_GLM;        % observation function


% Build priors structure
priors.muPhi = zeros(1,1);         % prior mean on observation params
priors.SigmaPhi = 1e0*eye(1); % prior covariance on observation params
priors.a_sigma = sigma;             % Jeffrey's prior
priors.b_sigma = 1e0;             % Jeffrey's prior
options.priors = priors;        % include priors in options structure
options.Laplace = 1;
options.updateHP = 0;
dim.n_phi = 1;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states

% fill in priors
[options,u,dim] = VBA_check(zeros(1e2,1),[],[],g_fname,dim,options);
options.verbose = 0;
options.DisplayWin = 0;
priors = options.priors;
priors2 = priors;

gridv = 2.^[-14:2:52];
N = 16;

inG.X = ones(1e2,1);
phi = 1e0;
options.inG = inG;

for i=1:N

    i
    
    % simulate data under full model
    [gx] = feval(g_fname,[],phi,[],inG);
    y = gx + sqrt(sigma.^-1)*randn(size(gx));

    
    % Invert full model
    options.priors = priors;
    [p1,o1] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
    F(i) = o1.F;
    
    % Invert reduced model
    priors2.SigmaPhi(1,1) = 0;
%     F0(i) = VB_SavageDickey(p1,priors,o1.F,dim,priors2);
    options.priors = priors2;
    [p0,o0] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
    F0(i) = o0.F;
    
    % Use Savage-Dickey ratio to obtain full posteriors
    for j=1:length(gridv)
        priors2.SigmaPhi(1,1) = gridv(j);
        options.priors = priors2;
        [p2,o2] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
        F1(i,j) = o2.F;
        dF(i,j) = F1(i,j) - F0(i);
        mu(i,j) = p2.muPhi;
        v(i,j) = sqrt(p2.SigmaPhi);
        [F11(i,j),p2] = VB_SavageDickey(p1,priors,o1.F,dim,priors2);
        dF1(i,j) = F11(i,j) - F0(i);
    end
    
    
    
end


bf = 1./(1+exp(-dF));
bf1 = 1./(1+exp(-dF1));

figure,errorbar(gridv,mean(bf,1),std(bf,[],1))
hold on
plot(gridv,mean(bf,1),'.')
errorbar(gridv,mean(bf1,1),std(bf1,[],1),'g')
plot(gridv,mean(bf1,1),'g.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('1-P(phi=0|y,m)')
figure,errorbar(gridv,mean(mu,1),std(mu,[],1))
hold on
plot(gridv,mean(mu,1),'.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('E[phi|y,m]')

figure,errorbar(gridv,mean(v,1),std(v,[],1))
hold on
plot(gridv,mean(v,1),'.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('V[phi|y,m]')


getSubplots








