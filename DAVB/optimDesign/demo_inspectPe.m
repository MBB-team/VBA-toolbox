% Demo: evaluation of design optimization for DCM comparison.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

close all
clear variables



%-----------------------------------------------------------
%-------------- common DCM model specification --------------------

%--------- Basic settings ----------------------------------
n_t = 2.5e1;                    % number of time samples
TR = 3e0;                       % sampling period (in sec)
microDT = 5e-2;                 % micro-time resolution (in sec)
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
homogeneous = 1;
reduced_f = 1;
lin = 1;
nreg = 2;
nu = 2;

%--------- model 1 -----------------------------------------
% invariant effective connectivity
A = [0 1
    1 0];
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0
    0 0];
% input-state coupling
C = [1 0
    0 1];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
% Build priors and options/dim structures for model inversion
[options{1},dim{1}] = getOptions4dcm(...
    A,...
    B,...
    C,...
    D,...
    TR,...
    microDT,...
    n_t,...
    homogeneous,...
    reduced_f,...
    lin);


%--------- model 2 -----------------------------------------
% invariant effective connectivity
A = [0 1
    1 0];
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0
    0 0];
% input-state coupling
C = [0 1
    1 0];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
% Build priors and options/dim structures for model inversion
[options{2},dim{2}] = getOptions4dcm(...
    A,...
    B,...
    C,...
    D,...
    TR,...
    microDT,...
    n_t,...
    homogeneous,...
    reduced_f,...
    lin);

%--------- design 1 ----------------------------------------
u{1} = zeros(2,n_t);
u{1}(1,1:5) = 1;
% u{1}(2,1:5) = 1;
% u{1}(1,16:20) = 1;
u{1}(2,16:20) = 1;

%--------- design 2 ----------------------------------------
u{2} = zeros(2,n_t);
u{2}(1,1:5) = 1;
% u{2}(2,1:5) = 1;
u{2}(1,16:20) = 1;
u{2}(2,16:20) = 1;


%--------- design 3 ----------------------------------------
u{3} = zeros(2,n_t);
u{3}(1,1:5) = 1;
u{3}(2,1:5) = 1;
u{3}(1,16:20) = 1;
u{3}(2,17:21) = 1;

% simulation parameters

nnp = 16;
pa = zeros(nnp,1);



% pa = {exp(-0.5);exp(-1.5)};
% pc = {exp(0.1);exp(0.1)};



sigma = {1e-1,5e-2};
indTrue = {1,2}; % 'true model' index


%-----------------------------------------------------------
%---------------- get design efficiency --------------------
nd = length(u); % number of designs
nm = length(dim); % number of models
ns = length(sigma); % number of noise precisions
np = length(pa); % number of simulation parameters sets
ni = length(indTrue); % 'true model' factor
NN = 16; % number of Monte-Carlo repetitions

Pe = zeros(nd,ns,np,ni);

if ~exist('i1234.mat','file')
    disp('initiate multi-session log file...')
    lo.i = 0;
    lo.i1 = 0;
    lo.i2 = 0;
    lo.Pe = Pe;
    save('i1234.mat','lo')
end

for i=1:nd %loop over experimental designs
    
    
    for i1=1:ns % noise precision factor
        
        disp(['Noise precision: ',num2str(sigma{i1})])
        
        for j=1:nm % loop over models
            % gets laplace approx to model j under design i
            options{j}.priors.b_sigma = 1./sigma{i1};
            [muy{i,j},Vy{i,j}] = getLaplace(...
                u{i},...
                f_fname,...
                g_fname,...
                dim{j},...
                options{j});
        end
        
        
        % computes Jensen-Shannon divergence of models under design i
        ps = ones(length(dim),1);
        ps = ps./sum(ps);
        [DJS(i,i1),b(i,i1)] = JensenShannon(muy(i,:),Vy(i,:),ps);
        
        
        
        for i2=1:np % prior mean to simulated parameter distance
            
            i2
            load('i1234.mat');
            
            if i>lo.i || i1>lo.i1 || i2>lo.i2
                lo.i = i;
                lo.i1 = i1;
                lo.i2 = i2;
                save('i1234.mat','lo')
                
                
                for i3=1:ni % 'true model' factor
                    
                    nerrors = 0;
                    
                    for i4=1:NN % Monte-Carlo repetitions
                        
                        
                        dim0 = dim{i3};
                        options0 = options{i3};
                        
                        stop = 0;
                        while ~stop
                            
                            %                         pa0 = pa(:,i2);
                            %                         pc0 = pa(:,i2);
                            
                            pa0 = sqrt(0.1)*randn(2,1);
                            pc0 = sqrt(0.1)*randn(2,1);
                            
                            [theta,phi] = getParamSImulDCM(...
                                homogeneous,...
                                dim0,...
                                options0,...
                                pa0,...
                                pc0);
                            
                            [y,x,x0,eta,e] = simulateNLSS(...
                                n_t,...
                                f_fname,...
                                g_fname,...
                                theta,...
                                phi,...
                                u{i},...
                                Inf,...
                                sigma{i1},...
                                options0);
                            
                            if ~isweird(y) && ~isweird(x)
                                stop = 1;
                            end
                            
                        end
                        
                        %                     % display time series of hidden states and observations
                        %                     displaySimulations(y,x,eta,e)
                        %                     pause
                        
                        for j=1:nm % loops over models
                            options{j}.DisplayWin = 0;
                            [posterior,out] = VBA_NLStateSpaceModel(...
                                y,...
                                u{i},...
                                f_fname,...
                                g_fname,...
                                dim{j},...
                                options{j});
                            F(j) = out.F;
                        end
                        win = find(F==max(F));
                        nerrors = nerrors+~isequal(win,i3);
                        
                    end
                    
                    Pe(i,i1,i2,i3) = nerrors./NN;
                    lo.Pe(i,i1,i2,i3) = Pe(i,i1,i2,i3);
                    save('i1234.mat','lo')
                    
                    
                end
                
                
            else
                
                disp(['skipping: i=',num2str(i),...
                    ' ; i1=',num2str(i1),...
                    ' ; i2=',num2str(i2)])
                
            end
            
        end
        
        
    end
    
    
    
    
end

% ep = mean(mean(Pe,4),3);
% ntot = NN*ni*np;
% a0 = ep.*ntot;
% b0 = ntot- a0;
% sp = sqrt(a0.*b0./(((a0+b0).^2).*(a0+b0+1)));
% 
% ub = b/2;
% lb = b.^2/4;
% 
% hf = figure('color',[1 1 1]);
% 
% for i=1:2
%     ha(i) = subplot(1,2,i,'parent',hf);
%     title(ha(i),['noise precision = ',num2str(sigma{i})])
%     set(ha(i),...
%         'nextplot','add',...
%         'xtick',1:nd,...
%         'xlim',[0.5,nd+0.5],...
%         'ylim',[-0.01 0.51],...
%         'ygrid','on')
%     
%     xt = cell(0);
%     for j = 1:size(b,1)
%         plot(ha(i),...
%             [j j],[ep(j,i)-sp(j,i),ep(j,i)+sp(j,i)],...
%             'color',0.5*[1 1 1],...
%             'linewidth',12)
%         plot(ha(i),...
%             j,ep(j,i),...
%             'k+',...
%             'linewidth',2)
%         if b(j,i) > 0
%             plot(ha(i),...
%                 [j j],[ub(j,i),lb(j,i)],...
%                 'r',...
%                 'linewidth',2)
%         else
%             plot(ha(i),...
%                 j,ep(j,i),...
%                 'ro',...
%                 'linewidth',2)
%         end
%         xt{j} = ['design ',num2str(j)];
%         
%     end
%     set(ha(i),'xticklabel',xt);
%     
% end












