% dummy susceptibility analysis
% this demo evaluates suscpetibility analyses in the following input-output
% scenario:
% o = P1*u1 + P2*(u1+u2) + u3*(P3+P4) + u4*P5 + P6.
% In other words:
% - u1 flows through two parallel branches: u1 --[P1 P2         ]--> o
% - u2 flows through one branch:            u2 --[   P2         ]--> o
% - u3 flows through two parallel branches: u3 --[      P3 P4   ]--> o
% - u4 flows through one branch:            u4 --[            P5]--> o
% This also means that:
% - parameters P1 and P2 are important for funnelling the u1->o
% relationship (and nothing else)
% - P2 is important for funnelling the u2->o relationship (and nothing
% else)
% - P3 and P4 are important for funnelling the u3->o relationship(and
% nothing else)
% - P5 is important for funnelling the u4->o relationship (and nothing
% else)
% - P6 is not involved in any input-output relationship
% Irrespective of the principles their stem from, any susceptibility
% analysis should "discover" this hidden structure.

clear variables
close all
clc

N = 1;


g_fname = @g_mixU;
inG.weights = [1;1;1;1;1;1];

dim.p       = 1;
dim.n       = 0;
dim.n_t     = 256;
dim.n_phi   = 6;
dim.n_theta = 0;

options.dim        = dim;
options.inG        = inG;
options.DisplayWin = 0;
options.verbose    = 0;

sigma              = 0.1;
phi                = [1;1;1;1;1;1];

for i=1:N
    disp(i)
    u = randn(4,dim.n_t);
    [y,x,x0,eta,e] = VBA_simulate (dim.n_t,[],g_fname,[],phi,u,[],sigma,options,zeros(dim.n,1));
    [posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
    % displayResults(posterior,out,y-e,[],[],[],phi,sigma,[])
    result(i) = VBA_susceptibility(posterior,out);
end

[result.susceptibility] ; susceptibility = [ans.norm] ;
[result.specificity]; specificity = [ans.norm] ;

xiPhi = cat(3,susceptibility.phi) ;
hf = figure('color',[1 1 1],'name','relative susceptibility');
ha = axes('parent',hf);
hi = imagesc(mean(xiPhi,3),'parent',ha);
colorbar('peer',ha)
colormap(flipud(colormap('bone')))
xlabel(ha,'parameters')
ylabel(ha,'inputs')
set(ha,'xtick',1:size(xiPhi,2),'ytick',1:size(xiPhi,1))
set(gca,'CLim',[0 1]) ;

zetaPhi = cat(3,specificity.phi) ;

hf = figure('color',[1 1 1],'name','relative specificity');
ha = axes('parent',hf);
hi = imagesc(mean(zetaPhi,3),'parent',ha);
colorbar('peer',ha)
colormap(flipud(colormap('bone')))
xlabel(ha,'parameters')
ylabel(ha,'inputs')
set(ha,'xtick',1:size(zetaPhi,2),'ytick',1:size(zetaPhi,1))
set(gca,'CLim',[0 1]) ;




%% interactions
figure
interaction = [result.interaction];



subplot(1,4,1)
inter_norm = [interaction.norm] ;
inter      = cat(3,inter_norm.phi) ;

imagesc(mean(inter,3));
colorbar()
colormap(flipud(colormap('bone')))
xlabel('parameters')
ylabel('inputs')
set(gca,'xtick',1:size(inter,2),'ytick',1:size(inter,1))
set(gca,'CLim',[0 1]) ;



subplot(1,4,2)
inter_norm = [interaction.normU] ;
inter      = cat(3,inter_norm.phi) ;

imagesc(mean(inter,3));
colorbar()
colormap(flipud(colormap('bone')))
xlabel('parameters')
ylabel('inputs')
set(gca,'xtick',1:size(inter,2),'ytick',1:size(inter,1))
set(gca,'CLim',[0 1]) ;



subplot(1,4,3)
inter_norm = [interaction.normP] ;
inter      = cat(3,inter_norm.phi) ;

imagesc(mean(inter,3));
colorbar()
colormap(flipud(colormap('bone')))
xlabel('parameters')
ylabel('inputs')
set(gca,'xtick',1:size(inter,2),'ytick',1:size(inter,1))
set(gca,'CLim',[0 1]) ;



subplot(1,4,4)
inter_norm = [interaction.normC] ;
inter      = cat(3,inter_norm.phi) ;

imagesc(mean(inter,3));
colorbar()
colormap(flipud(colormap('bone')))
xlabel('parameters')
ylabel('inputs')
set(gca,'xtick',1:size(inter,2),'ytick',1:size(inter,1))
set(gca,'CLim',[0 1]) ;

