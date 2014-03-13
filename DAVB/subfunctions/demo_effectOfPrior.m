% Demo for the multiple comparison problem

close all
clear variables

% Choose basic settings for simulations
sigma = 1e0;            % precision 
g_fname = @g_GLM;        % observation function
ny = 1e2;

% Build priors structure
priors.muPhi = 0;        % prior mean on observation params
priors.SigmaPhi = 1e0;     % prior covariance on observation params
priors.a_sigma = 1e0;             % Jeffrey's prior
priors.b_sigma = 1e0;             % Jeffrey's prior
options.priors = priors;        % include priors in options structure
dim.n_phi = 1;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states

% fill in priors
[options,u,dim] = VBA_check(zeros(ny,1),[],[],g_fname,dim,options);
options.verbose = 0;
options.DisplayWin = 0;
priors = options.priors;

gridv = 10.^[-6:.25:6];
N = 16;

inG.X = ones(ny,1);
options.inG = inG;

for i=1:N

    i
    
    % simulate data under full model...
    phi = 1e0;
    [gx] = feval(g_fname,[],phi,[],inG);
    y1 = gx + sqrt(sigma.^-1)*randn(size(gx));
    % ... and under the null
    phi = 0e0;
    [gx] = feval(g_fname,[],phi,[],inG);
    y0 = gx + sqrt(sigma.^-1)*randn(size(gx));
    
    % Derive Bayes' factor given both types of data
    for j=1:length(gridv)
        
        priors.SigmaPhi(1,1) = gridv(j);
        options.priors = priors;
        
        [p1,o1] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
        mu1(i,j) = p1.muPhi;
        v1(i,j) = sqrt(p1.SigmaPhi);
        
        [p0,o0] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,options);
        mu0(i,j) = p0.muPhi;
        v0(i,j) = sqrt(p0.SigmaPhi);
        
        priors.SigmaPhi(1,1) = 0;
        options.priors = priors;
        [p1,o10] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
        [p0,o00] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,options);
        
%         [F1(i,j)] = VBA_SavageDickey(p1,options.priors,o1.F,dim,priors);
        dF1(i,j) = o1.F - o10.F;%F1(i,j);
        
%         [F0(i,j)] = VBA_SavageDickey(p0,options.priors,o0.F,dim,priors);
        dF0(i,j) = o0.F - o00.F;%F0(i,j);
        
    end
    
    
    
end


bf1 = 1./(1+exp(-dF1));
bf0 = 1./(1+exp(-dF0));

figure
errorbar(gridv,mean(bf1,1),std(bf1,[],1))
hold on
errorbar(gridv,mean(bf0,1),std(bf0,[],1),'g')
plot(gridv,mean(bf1,1),'.')
plot(gridv,mean(bf0,1),'g.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('1-P(phi=0|y,m)')
legend('phi=1','phi=0')



figure
errorbar(gridv,mean(mu1,1),std(mu1,[],1))
hold on
errorbar(gridv,mean(mu0,1),std(mu0,[],1),'g')
plot(gridv,mean(mu1,1),'.')
plot(gridv,mean(mu0,1),'g.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('E[phi|y,m]')
legend('phi=1','phi=0')

figure
errorbar(gridv,mean(v1,1),std(v1,[],1))
hold on
errorbar(gridv,mean(v0,1),std(v0,[],1),'g')
plot(gridv,mean(v1,1),'.')
plot(gridv,mean(v0,1),'g.')
box off
axis tight
set(gca,'xscale','log')
xlabel('prior variance')
ylabel('V[phi|y,m]')
legend('phi=1','phi=0')

getSubplots








