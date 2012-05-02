% Demo: design optimization for DCM of fMRI data.


close all
clear variables



%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--------- Basic settings ----------------------------------
n_t = 2.5e1;                    % number of time samples
TR = 3e0;                       % sampling period (in sec)
microDT = 5e-2;                 % micro-time resolution (in sec)
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;


%--------- model 1 -----------------------------------------
% invariant effective connectivity
A = [0 1
     1 0];
nreg = size(A,1);
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
[options{1},dim{1}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t);


%--------- model 2 -----------------------------------------
% invariant effective connectivity
A = [0 1
     1 0];
nreg = size(A,1);
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
[options{2},dim{2}] = getOptions4dcm(A,B,C,D,TR,microDT,n_t);


%--------- design 1 ----------------------------------------
u{1} = zeros(2,n_t);
u{1}(1,1:5) = 1;
u{1}(2,1:5) = 1;
u{1}(1,16:20) = 1;
u{1}(2,16:20) = 1;

%--------- design 2 ----------------------------------------
u{2} = zeros(2,n_t);
u{2}(1,1:5) = 1;
% u{2}(2,1:5) = 1;
% u{2}(1,16:20) = 1;
u{2}(2,16:20) = 1;


%--------- design 3 ----------------------------------------
u{3} = zeros(2,n_t);
u{3}(1,1:5) = 1;
% u{1}(2,1:5) = 1;
u{3}(1,16:20) = 1;
u{3}(2,16:20) = 1;

%-----------------------------------------------------------
%---------------- get design efficiency --------------------
nd = length(u);
nm = length(dim);
for i=1:nd %loops over experimental designs
    
    fprintf(1,'\n')
    fprintf(1,['design ',num2str(i),' / '])
    
%     fprintf(1,['design ',num2str(i),'/',num2str(nd),' , model       '])
%     
%     for j=1:nm % loops over models
%         % gets laplace approx to model j under design i
%         fprintf(1,'\b\b\b\b\b\b')
%         fprintf(1,[num2str(j),'/',num2str(nm),'...'])
%         [muy{i,j},Vy{i,j}] = getLaplace(...
%             u{i},...
%             f_fname,...
%             g_fname,...
%             dim{j},...
%             options{j},...
%             1);
%     end
%     fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
%     fprintf(1,[' : OK.'])
%     
%     
%     % computes Jensen-Shannon divergence of models under design i
%     ps = ones(length(dim),1);
%     ps = ps./sum(ps);
%     [DJS(i),b(i)] = JensenShannon(muy(i,:),Vy(i,:),ps);
    
    [DJS(i),out] = designEfficiency(f_fname,g_fname,dim,options,u{i},'models');
    out
    b(i) = out.b;
    for j=1:nm % loops over models
        muy{i,j} = out.muy{j};
        Vy{i,j} = out.Vy{j};
    end
    
end
fprintf(1,'\n\n')

[h] = displayOptimDesign(muy,Vy,u,dim,DJS,b);


