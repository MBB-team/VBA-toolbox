% Demo for the multiple comparison problem

close all
clear variables

% Choose basic settings for simulations
sigma = 1e0;            % precision 
g_fname = @g_GLM;        % observation function
ny = 1e2;

% Build options/dim structure
options.verbose = 0;
options.DisplayWin = 0;
dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n = 0;                        % nb of hidden states
dim.n_t = 1;
dim.p = ny*2;

% form 2x2 model space using prior structures
p0.muPhi = zeros(2,1);        % prior mean on observation params
p0.SigmaPhi = 1e0*eye(2);     % prior covariance on observation params
p0.a_sigma = 1e0;             % Jeffrey's prior
p0.b_sigma = 1e0;             % Jeffrey's prior
for i=1:2
    for j=1:2
        priors{i,j} = p0;
        priors{i,j}.SigmaPhi(1,1) = 2-i;
        priors{i,j}.SigmaPhi(2,2) = 2-j;
    end
end

N = 16;

gridp = 0:5e-2:1;
X = [ones(ny,1);zeros(ny,1)];
X = [X,1-X];

% full model
phi1 = [1;1];
% null model
phi0 = [0;1];

for imcmc=1:N

    imcmc
    
  
    for k=1:length(gridp)
        
        % sample GLM regressors
        A = [[1;0],[gridp(k);1-gridp(k)]];
        inG.X = X*A;
        
        % simulate data under full model...
        [gx] = feval(g_fname,[],phi1,[],inG);
        e = sqrt(sigma.^-1)*randn(size(gx));
        e(1:ny) = e(1:ny)-mean(e(1:ny));
        e(ny+1:2*ny) = e(ny+1:2*ny)-mean(e(ny+1:2*ny));
        y1 = gx + e;
        % ... and under the null
        [gx] = feval(g_fname,[],phi0,[],inG);
        y0 = gx + e;
        
        % invert full model on both data sets
        options.inG = inG;
        options.priors = priors{1,1};
        [p1,o1] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
        [p0,o0] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,options);
        
        % Use Savage-Dickey ratio for both types of data
        for i=1:2
            for j=1:2
                if i==1 && j==1
                    F1(i,j) = o1.F;
                    F0(i,j) = o0.F;
                else
                    options.priors = priors{i,j};
                    [p1,o10] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,options);
                    [p0,o00] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,options);
                    F1(i,j) = o10.F;
                    F0(i,j) = o00.F;
                    
%                     [F1(i,j)] = VBA_SavageDickey(p1,o1.options.priors,o1.F,dim,priors{i,j});
%                     [F0(i,j)] = VBA_SavageDickey(p0,o0.options.priors,o0.F,dim,priors{i,j});
                end
            end
        end
        P1(:,:,k,imcmc) = exp(F1);
        P1(:,:,k,imcmc) = P1(:,:,k,imcmc)./sum(sum(P1(:,:,k,imcmc)));
        P0(:,:,k,imcmc) = exp(F0);
        P0(:,:,k,imcmc) = P0(:,:,k,imcmc)./sum(sum(P0(:,:,k,imcmc)));
        
    end
        
end


hf = figure;
for i=1:4
    ha(i) = subplot(2,2,i,'parent',hf,'nextplot','add');
end
for k=1:length(gridp)
    
    m = mean(P1(1,1,k,:),4);
    sd = std(P1(1,1,k,:),[],4);
    errorbar(ha(1),gridp(k),m,sd)
    plot(ha(1),gridp(k),m,'.')
    plot(ha(1),[gridp(1),gridp(end)],[0.95,0.95],'r')
    plot(ha(1),[gridp(1),gridp(end)],[0.05,0.05],'r')
    title(ha(1),'full model space')
    
    P1m = P1(1,1,k,:) + P1(1,2,k,:); % proba for Phi(1)~=0
    P2m = P1(1,1,k,:) + P1(2,1,k,:); % proba for Phi(2)~=0
    pp = P1m.*P2m; % proba for Phi(1) and Phi(2) ~=0
    m = mean(pp,4);
    sd = std(pp,[],4);
    errorbar(ha(2),gridp(k),m,sd)
    plot(ha(2),gridp(k),m,'.')
    plot(ha(2),[gridp(1),gridp(end)],[0.95,0.95],'r')
    plot(ha(2),[gridp(1),gridp(end)],[0.05,0.05],'r')
    title(ha(2),'joint inference')
    
    m = mean(P1m,4);
    sd = std(P1m,[],4);
    errorbar(ha(3),gridp(k),m,sd)
    plot(ha(3),gridp(k),m,'.')
    plot(ha(3),[gridp(1),gridp(end)],[0.95,0.95],'r')
    plot(ha(3),[gridp(1),gridp(end)],[0.05,0.05],'r')
    title(ha(3),'testing for Phi(1)')
    
    m = mean(P2m,4);
    sd = std(P2m,[],4);
    errorbar(ha(4),gridp(k),m,sd)
    plot(ha(4),gridp(k),m,'.')
    plot(ha(4),[gridp(1),gridp(end)],[0.95,0.95],'r')
    plot(ha(4),[gridp(1),gridp(end)],[0.05,0.05],'r')
    title(ha(4),'testing for Phi(2)')
    
    
    
    % under the null now:
    
    m = mean(P0(1,1,k,:),4);
    sd = std(P0(1,1,k,:),[],4);
    errorbar(ha(1),gridp(k),m,sd,'g')
    plot(ha(1),gridp(k),m,'g.')
    
    P1m = P0(1,1,k,:) + P0(1,2,k,:); % proba for Phi(1)~=0
    P2m = P0(1,1,k,:) + P0(2,1,k,:); % proba for Phi(2)~=0
    pp = P1m.*P2m; % proba for Phi(1) and Phi(2) ~=0
    m = mean(pp,4);
    sd = std(pp,[],4);
    errorbar(ha(2),gridp(k),m,sd,'g')
    plot(ha(2),gridp(k),m,'g.')
    
    m = mean(P1m,4);
    sd = std(P1m,[],4);
    errorbar(ha(3),gridp(k),m,sd,'g')
    plot(ha(3),gridp(k),m,'g.')
    
    m = mean(P2m,4);
    sd = std(P2m,[],4);
    errorbar(ha(4),gridp(k),m,sd,'g')
    plot(ha(4),gridp(k),m,'g.')
    
    
end

getSubplots








