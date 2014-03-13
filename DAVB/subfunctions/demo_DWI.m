% Demo for Gaussian convolution observation model.

close all
clear variables

%---- simulate noisy gbf ----%

% Choose basic settings for simulations
sigma = 1e5;            % precision
n = 10;

xg = 1:n;
yg = 1:n;
[X,Y] = meshgrid(xg,yg);
inG.P = [ones(n^2,1),vec(X),vec(Y),vec(X).*vec(Y),vec(X).^2,vec(Y).^2];
inG.ana = 1;
phi = zeros(6*9,1);
ind_diag = find(vec(eye(3)));
id = [];
for i=1:9
    ind = (i-1)*6+1:i*6;
    if ismember(i,ind_diag)
        sc = 10;
        mu = -0.5;
        id = [id;ind(:)];
    else
        sc = 1;
        mu = -0.5;
    end
    
    Phi(:,i) = sc*(rand(6,1)+mu);
    phi(ind) = Phi(:,i);
end




% I = speye(n,n);
% E = sparse(2:n,1:n-1,1,n,n);
% D = E+E'-2*I;
% A = kron(D,I)+kron(I,D);
% [u,s,v] = svd(full(A));
% L = u*diag(1./diag(s))*v';
% L2 = u*diag(1./(diag(s).^2))*v';
% 
% for i=1:9
%     tmp = randn(n,n);
%     Phi(:,i) = L*vec(tmp);
% end
% phi = vec(Phi);
   
g_fname = @g_DWI; % observation function

% Build simulated observations 
[gx] = feval(g_fname,[],phi,[],inG);
y = gx + sqrt(sigma.^-1)*randn(size(gx));
% display time series of hidden states and observations
figure,
plot(y','ro')
hold on;
plot(gx')
% disp('--paused--')
% pause

hf = figure;
ha = axes('parent',hf);
for i=1:6
    ind = (i-1)*n^2+1:i*n^2;
    Iy(:,:,i) = reshape(gx(ind),n,n);
    imagesc(Iy(:,:,i),'parent',ha)
    pause
end
% 
% for i=1:9
%     ind = (i-1)*6+1:i*6;
%     vecA = inG.P*phi(ind);
%     Iphi(:,:,i) = reshape(vecA,n,n);
%     imagesc(Iphi(:,:,i),'parent',ha)
%     pause
% end


%---- Invert gbf on simulated data ----%

% Build priors structure
% priors.muPhi = zeros(9*n.^2,1);         % prior mean on observation params
% priors.SigmaPhi = 1e0*kron(eye(9),L2);  % prior covariance on observation params
priors.muPhi = 1e-3*ones(6*9,1);
% priors.muPhi(ind_diag) = 6;
priors.SigmaPhi = 1e-4*eye(6*9);
priors.SigmaPhi(ind_diag,id) = 1e-2;
priors.a_sigma = 1e5;                   % Jeffrey's prior
priors.b_sigma = 1e0;                     % Jeffrey's prior
% Build options structure
options.priors = priors;        % include priors in options structure
options.GnFigs = 0;             % disable annoying figures
options.verbose = 1;
options.inG = inG;
% dim.n_phi = n.^2*9;                  % nb of observation parameters
dim.n_phi = 6*9;
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],phi,[],sigma)

hf = figure;
for i=1:4
    ha(i) = subplot(2,2,i,'parent',hf);
end
for i=1:9
    ind = (i-1)*6+1:i*6;
    vecA = inG.P*phi(ind);
    Iphi(:,:,i) = reshape(vecA,n,n);
    vecAhat = inG.P*posterior.muPhi(ind);
    Iphihat(:,:,i) = reshape(vecAhat,n,n);
    
    imagesc(Iphi(:,:,i),'parent',ha(1))
    title(ha(1),'simulated field')
    colorbar('peer',ha(1))
    
    imagesc(Iphihat(:,:,i),'parent',ha(2))
    title(ha(2),'estimated field')
    colorbar('peer',ha(2))
    
    imagesc(Iphi(:,:,i)-Iphihat(:,:,i),'parent',ha(3))
    title(ha(3),'estimated-simulated field')
    colorbar('peer',ha(3))
    
    plot(vec(Iphi(:,:,i)),vec(Iphihat(:,:,i)),'.','parent',ha(4))
    title(ha(4),'estimated vs simulated field')
    
    pause
end



