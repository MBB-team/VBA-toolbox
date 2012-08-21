% This demo evaluates the ability of VB to accept the null using MCMC
% simulation series of data generated under H0 vs under H1


clear all
close all


beta = 1;
s2 = 1e-1;
n = 10;

N = 1e2;

options.inG.X = ones(n,1);
options.DisplayWin = 0;
options.verbose =0;
options.priors.muPhi = 0;
options.priors.SigmaPhi = 1;
options.Laplace = 1;
dim.n_phi = 1;
dim.n = 0;
dim.n_theta = 0;
g_fname = @g_GLM;
o1 = options;
o0 = options;
o0.priors.SigmaPhi = 0;

LBF = zeros(N,2);

for i=1:N
    
    i
    
    y0 = 0 + sqrt(s2).*randn(n,1);
    [p1,ou1] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,o1);
    [p0,ou0] = VBA_NLStateSpaceModel(y0,[],[],g_fname,dim,o0);
    LBF(i,1) = ou1.F - ou0.F;

    y1 = beta + sqrt(s2).*randn(n,1);
    [p1,ou1] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,o1);
    [p0,ou0] = VBA_NLStateSpaceModel(y1,[],[],g_fname,dim,o0);
    LBF(i,2) = ou1.F - ou0.F;
    
end

[ny,nx] = hist(LBF,30);
figure('color',[1 1 1 ])
bar(nx,ny)
hold on
plot([-3,-3],get(gca,'ylim'),'k')
plot([3,3],get(gca,'ylim'),'k')
legend({'y~p(y|H0)','y~p(y|H1)'})
xlabel('LBF = log p(y|H1) - log p(y|H0)')



