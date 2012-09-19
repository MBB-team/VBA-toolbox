% Demo for Cramer-Rao bound with non-identifiable parameters

% close all
clear variables

%---- simulate noisy data ----%

ny = 16;

gridp = 2.^[0:8];
N = 16;


for i=1:N
    
    i
    for j=1:length(gridp)
        
        
        nphi = gridp(j);
        
        % Choose basic settings for simulations
        sigma = 1e0;            % precision
        phi = ones(nphi,1);       % observation parameters
        g_fname = @g_GLM;        % observation function
        inG.X = randn(ny+1,nphi);
        % Build simulated observations
        [gx] = feval(g_fname,[],phi,[],inG);
        y = gx + sqrt(sigma.^-1)*randn(size(gx));
        % display time series of hidden states and observations
%         figure,
%         plot(y','ro')
%         hold on;
%         plot(gx')
        % disp('--paused--')
        % pause
        
        
        %---- Invert model on train data ----%
        
        ytrain = y(1:ny);
        options.inG.X = inG.X(1:ny,:);
        
        % Build priors structure
        priors.muPhi = zeros(nphi,1);         % prior mean on observation params
        priors.SigmaPhi = 1e1*eye(nphi); % prior covariance on observation params
        priors.a_sigma = 1e0;             % Jeffrey's prior
        priors.b_sigma = 1e0;             % Jeffrey's prior
        % Build options structure
        options.priors = priors;        % include priors in options structure
        options.verbose = 0;
        options.DisplayWin = 0;
        dim.n_phi = nphi;                  % nb of observation parameters
        dim.n_theta = 0;                % nb of evolution parameters
        dim.n=0;                        % nb of hidden states
        dim.p = ny;
        dim.n_t = 1;
        
        % Call inversion routine
        [posterior,out] = VBA_NLStateSpaceModel(ytrain,[],[],g_fname,dim,options);
        
        
        %--- Make prediction on test data ---%
        ytest = y(ny+1);
        opt.priors = posterior;
        opt.priors.iQy{1} = 1;
        opt.inG.X = inG.X(ny+1,:);
        
        dim.p = 1;
        dim.n_t = 1;
        [muy,Vy] = VBA_getLaplace([],[],g_fname,dim,opt);
        m = reshape(muy,dim.p,[]);
        v = reshape(diag(Vy),dim.p,[]);
        
        
%         hf = figure(...
%             'name','Predictive density: hidden states',...
%             'color',[1 1 1],...
%             'menubar','none');
%         hs = axes('parent',hf,'nextplot','add');
%         [haf,hf] = plotUncertainTimeSeries(m,v,[],hs);
%         plot(hs,ytest,'go')
%         pause
%         close(gcf)
        
        % store prediction error and expected prediction error
        
        PE(i,j) = (m-ytest).^2;
        EPE(i,j) = v;
        
        
    end
end


figure

errorbar(mean(PE,1),std(PE,[],1)./sqrt(N))
hold on
plot(mean(EPE),'r')
xlabel('(log) number of unknown variables')
legend({'observed squared (prediction) error','expected squared error'})
