% this demo exemplifies the use of mixed-effects analysis in VBA

clear all
close all
clc
dbstop if error

ns = 8; % #subjects
dim.n_phi = 2;
dim.n = 1;
dim.n_theta = 2;
dim.p = 1; % data dim (within-subject)
dim.n_t = 16; % trials (within-subject)

% simulate MFX
y = cell(ns,1);
nu = ones(5,1); % population mean
alpha = 100; % population precision
sigmas = ones(ns,1); % within-subject residual variance
g_fname = @g_vgo;
f_fname = @f_vgo;
for i=1:ns
    % draw within-subject effects from population distribution
    params(:,i) = nu + randn(5,1)./sqrt(alpha);
    theta(:,i) = params(1:2,i);
    phi(:,i) = params(3:4,i);
    x0(:,i) = params(5,i);
    u{i} = randn(1,dim.n_t);
    options{i}.DisplayWin = 0;
    options{i}.verbose = 0;
    options{i}.dim = dim;
    options{i}.binomial = 0;
    [y{i}] = simulateNLSS(dim.n_t,f_fname,g_fname,theta(:,i),phi(:,i),u{i},Inf,Inf,options{i},x0(:,i));
end


% priors_group.QPhi = 0.*eye(dim.n_phi);
% priors_group.QTheta = 0.*eye(dim.n_theta);
% priors_group.QX0 = 0.*eye(dim.n);
% priors_group.QPhi(2,2) = 0; % ffx
% priors_group.SigmaPhi = eye(dim.n_phi);
% priors_group.SigmaPhi(2,2) = 0; % fix population mean to 0
[p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_fname,g_fname,dim,options);%,priors_group);

% extract within-subject parameter estimates
for i=1:ns
    % with MFX-type priors 
    Theta(:,i,1) = p_sub{i}.muTheta;
    Phi(:,i,1) = p_sub{i}.muPhi;
    X0(:,i,1) = p_sub{i}.muX0;
    % without MFX-type priors
    Theta(:,i,2) = o_group.initVBA.p_sub{i}.muTheta;
    Phi(:,i,2) = o_group.initVBA.p_sub{i}.muPhi;
    X0(:,i,2) = o_group.initVBA.p_sub{i}.muX0;
end

hf = figure('name','btw-subject variability','color',[1 1 1]);
col = getColors(max([dim.n;dim.n_theta;dim.n_phi]));
tistr = {' (with MFX priors)',' (without MFX priors)'};
for j=1:2
    ha = subplot(2,3,(j-1)*3+1,'parent',hf,'nextplot','add');
    plot(ha,theta',Theta(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n_theta
        hp = plot(ha,theta(i,:),Theta(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
        if oa.pv<0.05
            leg{i} = ['dim #',num2str(i),': *'];
        else
            leg{i} = ['dim #',num2str(i)];
        end
    end
    axis(ha,'tight')
    title(ha,['evolution params',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
    ha = subplot(2,3,(j-1)*3+2,'parent',hf,'nextplot','add');
    plot(ha,phi',Phi(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n_phi
        hp = plot(ha,phi(i,:),Phi(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
        if oa.pv<0.05
            leg{i} = ['dim #',num2str(i),': *'];
        else
            leg{i} = ['dim #',num2str(i)];
        end
    end
    title(ha,['observation params',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
    ha = subplot(2,3,(j-1)*3+3,'parent',hf,'nextplot','add');
    plot(ha,x0',X0(:,:,j)','.')
    leg = cell(0);
    for i=1:dim.n
        hp = plot(ha,x0(i,:),X0(i,:,j),'.','color',col(i,:));
        [oa] = addBestLinearPredictor(hp,0);
        if oa.pv<0.05
            leg{i} = ['dim #',num2str(i),': *'];
        else
            leg{i} = ['dim #',num2str(i)];
        end
    end
    title(ha,['initial conditions',tistr{j}])
    xlabel(ha,'simulated')
    ylabel(ha,'estimated')
    legend(ha,leg)
end

