% Demo: evaluation of design optimization for DCM comparison: TMS?
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

close all
clear variables



%-----------------------------------------------------------
%-------------- common DCM model specification --------------------

%--------- Basic settings ----------------------------------
n_t =3e2;                    % number of time samples
TR = 1e0;                       % sampling period (in sec)
microDT = 1e-1;                 % micro-time resolution (in sec)
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
homogeneous = 1;
reduced_f = 0;
lin = 0;
nreg = 2;
nu = 2+4;
fbk = 1; % include models with feedback connections

%--------- model 1 -----------------------------------------
% invariant effective connectivity
A = [0 0
    0 0];
% modulatory effects
B{1} = [0 0
    0 0];
B{2} = [0 0
    0 0];
B{3} = [0 0
    0 0]; % TMS on-line on region 1
B{4} = [0 0
    0 0]; % TMS on-line on region 2
B{5} = [1 0
    0 0]; % TMS off-line on region 1
B{6} = [0 0
    0 1]; % TMS off-line on region 2
% input-state coupling
C = [1 0 1 0 0 0
    0 1 0 1 0 0];
% gating (nonlinear) effects
D{1} = [0 0
    0 1];
D{2} = [0 0
    0 0];
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
options{1}.priors.muTheta(1:options{1}.inF.indself) = 1e-1;
options{1}.priors.muTheta(1:options{1}.inF.indD{1}) = 5e-1;
tms_on = options{1}.inF.indC(3:4);
tms_off = cell2mat(options{1}.inF.indB(5:6));
options{1}.priors.muTheta(tms_on) = 1e-1;
options{1}.priors.muTheta(tms_off) = -1e1;
% options{1}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;


%--------- model 2 -----------------------------------------
% invariant effective connectivity
A = [0 0
    1 0];
% modulatory effects
B{1} = [0 0
    0 0];
B{2} = [0 0
    1 0];
B{3} = [0 0
    0 0]; % TMS on-line on region 1
B{4} = [0 0
    0 0]; % TMS on-line on region 2
B{5} = [1 0
    0 0]; % TMS off-line on region 1
B{6} = [0 0
    0 1]; % TMS off-line on region 2
% input-state coupling
C = [1 0 1 0 0 0
    0 0 0 1 0 0];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = [0 0
    0 0];
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
options{2}.priors.muTheta(1:options{2}.inF.indself) = 1e-1;
tms_on = options{2}.inF.indC(2:3);
tms_off = cell2mat(options{2}.inF.indB(5:6));
options{2}.priors.muTheta(tms_on) = 1e-1;
options{2}.priors.muTheta(tms_off) = -1e1;
% options{2}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;


%--------- model 3 -----------------------------------------
% invariant effective connectivity
A = [0 0
    1 0];
% modulatory effects
B{1} = [0 0
    0 0];
B{2} = [0 0
    0 1];
B{3} = [0 0
    0 0]; % TMS on-line on region 1
B{4} = [0 0
    0 0]; % TMS on-line on region 2
B{5} = [1 0
    0 0]; % TMS off-line on region 1
B{6} = [0 0
    0 1]; % TMS off-line on region 2
% input-state coupling
C = [1 0 1 0 0 0
    0 0 0 1 0 0];
% gating (nonlinear) effects
D{1} = [0 0
    0 0];
D{2} = [0 0
    0 0];
% Build priors and options/dim structures for model inversion
[options{3},dim{3}] = getOptions4dcm(...
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
options{3}.priors.muTheta(1:options{3}.inF.indself) = 1e-1;
options{3}.priors.muTheta(1:options{3}.inF.indB{2}) = 5e-1;
tms_on = options{3}.inF.indC(2:3);
tms_off = cell2mat(options{3}.inF.indB(5:6));
options{3}.priors.muTheta(tms_on) = 1e-1;
options{3}.priors.muTheta(tms_off) = -1e1;
% options{3}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;

if fbk
    %--------- model 4 -----------------------------------------
    % invariant effective connectivity
    A = [0 1
        0 0];
    % modulatory effects
    B{1} = [0 0
        0 0];
    B{2} = [0 0
        0 0];
    B{3} = [0 0
        0 0]; % TMS on-line on region 1
    B{4} = [0 0
        0 0]; % TMS on-line on region 2
    B{5} = [1 0
        0 0]; % TMS off-line on region 1
    B{6} = [0 0
        0 1]; % TMS off-line on region 2
    % input-state coupling
    C = [1 0 1 0 0 0
        0 1 0 1 0 0];
    % gating (nonlinear) effects
    D{1} = [0 0
        0 1];
    D{2} = [0 0
        0 0];
    % Build priors and options/dim structures for model inversion
    [options{4},dim{4}] = getOptions4dcm(...
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
    options{4}.priors.muTheta(1:options{4}.inF.indself) = 1e-1;
    options{4}.priors.muTheta(1:options{4}.inF.indD{1}) = 5e-1;
    tms_on = options{4}.inF.indC(3:4);
    tms_off = cell2mat(options{4}.inF.indB(5:6));
    options{4}.priors.muTheta(tms_on) = 1e-1;
    options{4}.priors.muTheta(tms_off) = -1e1;
    % options{1}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;
    
    
    %--------- model 5 -----------------------------------------
    % invariant effective connectivity
    A = [0 1
        1 0];
    % modulatory effects
    B{1} = [0 0
        0 0];
    B{2} = [0 0
        1 0];
    B{3} = [0 0
        0 0]; % TMS on-line on region 1
    B{4} = [0 0
        0 0]; % TMS on-line on region 2
    B{5} = [1 0
        0 0]; % TMS off-line on region 1
    B{6} = [0 0
        0 1]; % TMS off-line on region 2
    % input-state coupling
    C = [1 0 1 0 0 0
        0 0 0 1 0 0];
    % gating (nonlinear) effects
    D{1} = [0 0
        0 0];
    D{2} = [0 0
        0 0];
    % Build priors and options/dim structures for model inversion
    [options{5},dim{5}] = getOptions4dcm(...
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
    options{5}.priors.muTheta(1:options{5}.inF.indself) = 1e-1;
    tms_on = options{5}.inF.indC(2:3);
    tms_off = cell2mat(options{5}.inF.indB(5:6));
    options{5}.priors.muTheta(tms_on) = 1e-1;
    options{5}.priors.muTheta(tms_off) = -1e1;
    % options{2}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;
    
    
    %--------- model 6 -----------------------------------------
    % invariant effective connectivity
    A = [0 1
        1 0];
    % modulatory effects
    B{1} = [0 0
        0 0];
    B{2} = [0 0
        0 1];
    B{3} = [0 0
        0 0]; % TMS on-line on region 1
    B{4} = [0 0
        0 0]; % TMS on-line on region 2
    B{5} = [1 0
        0 0]; % TMS off-line on region 1
    B{6} = [0 0
        0 1]; % TMS off-line on region 2
    % input-state coupling
    C = [1 0 1 0 0 0
        0 0 0 1 0 0];
    % gating (nonlinear) effects
    D{1} = [0 0
        0 0];
    D{2} = [0 0
        0 0];
    % Build priors and options/dim structures for model inversion
    [options{6},dim{6}] = getOptions4dcm(...
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
    options{6}.priors.muTheta(1:options{6}.inF.indself) = 1e-1;
    options{6}.priors.muTheta(1:options{6}.inF.indB{2}) = 5e-1;
    tms_on = options{6}.inF.indC(2:3);
    tms_off = cell2mat(options{6}.inF.indB(5:6));
    options{6}.priors.muTheta(tms_on) = 1e-1;
    options{6}.priors.muTheta(tms_off) = -1e1;
    % options{3}.priors.SigmaTheta([tms_off;tms_on'],[tms_off;tms_on']) = 0;
    
    
end



N = 16; % # Monte-Carlo repetitions

for ii=1:N
    
    ii
    
    
    % design 1: TMS on-line on region 1
    u{1} = sample_u(3,16,0,1,ceil(n_t./TR),16,16);
    u{1} = [u{1};zeros(3,size(u{1},2))];
    
    % design 2: TMS on-line on region 2
    u{2} = sample_u(3,16,0,1,ceil(n_t./TR),16,16);
    u{2} = [u{2};zeros(3,size(u{2},2))];
    u{2}(4,:) = u{2}(3,:);
    u{2}(3,:) = 0;
    
    % design 3: TMS off-line on region 1
    u{3} = sample_u(2,16,0,1,ceil(n_t./TR),16,16);
    u{3} = [u{3};zeros(4,size(u{3},2))];
    u{3}(5,:) = 1;
    
    % design 4: TMS off-line on region 2
    u{4} = sample_u(2,16,0,1,ceil(n_t./TR),16,16);
    u{4} = [u{4};zeros(4,size(u{4},2))];
    u{4}(6,:) = 1;
    
    % design 5: no TMS
    u{5} = sample_u(2,16,0,1,ceil(n_t./TR),16,16);
    u{5} = [u{5};zeros(4,size(u{5},2))];
    
    
    
    
    %-----------------------------------------------------------
    %---------------- get design efficiency --------------------
    nd = length(u); % number of designs
    nm = length(dim); % number of models
    
    
    for i=1:nd %loop over experimental designs
        
        
        i
        
        for j=1:nm % loop over models
            % gets laplace approx to model j under design i
            %             options{j}.priors.b_sigma = 1./sigma{i1};
            %             options{j}.priors.muTheta(1:options{j}.inF.indself) = 1e-1;
            
            options{j}.priors.a_sigma = 5e-2;
            options{j}.priors.b_sigma = 1e0;
            
            
            [muy{i,j},Vy{i,j}] = getLaplace(...
                u{i},...
                f_fname,...
                g_fname,...
                dim{j},...
                options{j});
            
            %
            %             [y,x,x0,eta,e] = simulateNLSS(...
            %                 n_t,...
            %                 f_fname,...
            %                 g_fname,...
            %                 options{j}.priors.muTheta,...
            %                 options{j}.priors.muPhi,...
            %                 u{i},...
            %                 Inf,...
            %                 Inf,...
            %                 options{j});
            %
            %             displaySimulations(y,x,eta,e)
            %             set(gcf,'name',['model ',num2str(j),' ; design ',num2str(i)])
            %             subplot(3,2,6),imagesc(u{i})
            %             pause
            
            
        end
        
        
        % computes Jensen-Shannon divergence of models under design i
        ps = ones(length(dim),1);
        ps = ps./sum(ps);
        [DJS(i,ii),b(i,ii)] = JensenShannon(muy(i,:),Vy(i,:),ps);
        
        
    end
    
    
    
    
end







