function [out] = mediationAnalysis0(Y,X,M,options)
% performs a simple mediation analysis
% function [out] = mediationAnalysis0(Y,X,M)
% Let us consider the following regression equations:
% M = X*a + e_M
% Y = M*b + X*c + e_Y
% where X is the independent variable, M is the mediator variable and Y is
% the dependent variable.
% Testing whether the effect of X onto Y is mediated by M reduces to
% testing the indirect path ab=0, or equivalently: a=0 OR b=0. This is the
% statistical test performed by this function.
% This function first performs Baron & Kenny's 3-steps mediation analysis:
% 1- regress Y on X
% 2- regress M on X
% 3- regress Y on both X and M
% Then, it performs the Sobel test, i.e. it asks whether the relationship
% between X and Y is significantly reduced when including M. This reduction
% is due to the indirect pathway.
% NB: all effect sizes are reported in terms of adjusted percentage of
% variance explained!
% IN:
%   - Y: the dependent variable
%   - X: the independent variable
%   - M: the mediator variable
%   - options: an optional structure containing the following fields:
%       .alpha: significance level {0.05}
%       .verbose: verbose mode {1}
%       .display: display result in matlab window {1}
% OUT:
%   - out: output structure containing the following fields:
%       .b: 2X2 array of regression effects:
%           [   Y=bX  Y-cM=bX
%               Y=bM  Y-cX=bM  ]
%       .p: 2X2 array of corresponding p-values, i.e. P(b=0|H0)
%           [   P(X->Y|H0)  P(X->Y|M,H0)
%               P(M->Y|H0)  P(M->Y|X,H0)  ]
%       .sobel: structure containing the following fields:
%           .ab: the estimated effect of the indirect path [X->M->Y]
%           .sab: the approximate standard error on ab
%           .p: Sobel's p-value, i.e. P(ab=0|H0)
%       .montecarlo: structure containing the following fields:
%           .N: number of Monte-Carlo samples
%           .p: p-value, i.e. P(ab=0|H0)
%           .IC: confidence interval for ab (at significance level alpha)
%       .conj: structure containing the following fields:
%           .p: conjunctive-null p-value
%           .pa: p-value for a=0
%           .pb: p-value for b=0

n = size(Y,1);
if size(X,1)~= n || size(M,1)~= n
    disp('Error: input variables must have the same size!')
    out = [];
    return
end

try;alpha=options.alpha;catch;alpha=0.05;end
try;verbose=options.verbose;catch;verbose=1;end
try;display=options.display;catch;display=1;end


% regress Y on X
[p1,stat1,df1,a1] = GLM_contrast([ones(n,1),X],Y,[0;1],'F',0);

% regress M on X
[p2,stat2,df2,a2] = GLM_contrast([ones(n,1),X],M,[0;1],'F',0);

% regress Y on both X and M
[p3,stat3,df3,a3] = GLM_contrast([ones(n,1),X,M],Y,[0;1;0],'F',0);
[p4,stat4,df4,a4] = GLM_contrast([ones(n,1),X,M],Y,[0;0;1],'F',0);

% path analysis
de = a3.R2_a;
ie = a2.R2_a.*a4.R2_a;
te = a4.R2-a4.R2_a+ie;

% Sobel test
sobel.a = a2.b(2);
sobel.va = a2.vhat.*a2.iC(2,2);
sobel.b = a3.b(3);
sobel.vb = a3.vhat.*a3.iC(3,3);
sobel.ab = sobel.a.*sobel.b;
sobel.sab = sqrt(sobel.a.^2*sobel.vb + sobel.b^2*sobel.va);
sobel.p = 2*VBA_spm_Ncdf(-abs(sobel.ab./sobel.sab),0,1);

% Monte-Carlo test
montecarlo.N = 1e6; % # Monte-Carlo samples
A = sqrt(sobel.va).*randn(montecarlo.N,1);
B = sqrt(sobel.vb).*randn(montecarlo.N,1);
usAB = (abs(sobel.a)+A).*(abs(sobel.b)+B);
sAB = (sobel.a+A).*(sobel.b+B);
montecarlo.p = min([1,2*(1-length(find(usAB>0))./montecarlo.N)]);
montecarlo.IC = [VBA_quantile(sAB,alpha/2),VBA_quantile(sAB,1-alpha/2)];

% conjunctive test
conj.p = max([p2,p4]);
conj.pa = p2;
conj.pb = p4;

% gridz = 0:1e-3:10;
% [gridp] = normcdf(-gridz,0,1) + 1-normcdf(gridz,0,1);
% IC = [sobel.ab-gridz.*sobel.sab;sobel.ab+gridz.*sobel.sab];
% in = IC(1,:)<=0&IC(2,:)>=0;
% sobel.p = gridp(find(in==1,1));
% if isempty(sobel.p)
%     sobel.p = 0;
% end



% summary (Baron & Kenny)
if p2 > alpha
    os = 'no mediation';
else
    if p4 > alpha
        os = 'no mediation';
    else
        if sobel.p < alpha && p1 < alpha
            if p3 < alpha
                os = 'partial mediation';
            else
                os = 'full mediation';
            end
        else
            os = 'no conclusion';
        end
    end
end

% store results
out.p = [[p1;p2],[p3;p4]];
out.b = [[a1.b(2);a2.b(2)],VBA_vec(a3.b(2:3))];
out.sobel = sobel;
out.montecarlo = montecarlo;
out.conj = conj;
out.summary = os;


if verbose
    disp(' ')
    disp(['Date: ',datestr(clock)])
    disp(' ')
    disp('-- Regression results --')
    disp(['step 1) regress Y on X: R2=',num2str(round(1e3*a1.R2_a)/10),'% (p=',num2str(out.p(1,1)),')'])
    disp(['step 2) regress M on X: R2=',num2str(round(1e3*a2.R2_a)/10),'% (p=',num2str(out.p(2,1)),')'])
    disp(['step 3) regress Y on both X and M:'])
    disp(['   - X: R2[adj]=',num2str(round(1e3*a3.R2_a)/10),'% (p=',num2str(out.p(1,2)),')'])
    disp(['   - M: R2[adj]=',num2str(round(1e3*a4.R2_a)/10),'% (p=',num2str(out.p(2,2)),')'])
    disp(' ')
    disp('-- Path analysis --')
    disp(['Indirect effect X->M->Y: R2=',num2str(round(1e3*ie)/10),'% (Sobel test: p=',num2str(out.sobel.p),')'])
    disp(['Direct effect X->Y: R2=',num2str(round(1e3*de)/10),'% (p=',num2str(out.p(1,2)),')'])
    disp(['Total effect: R2=',num2str(round(1e3*te)/10),'%'])
    disp(['Summary: ',os])
    disp(' ')
    disp('-- Monte-Carlo sampling --')
    disp(['Indirect effect X->M->Y: P(ab=0|H0)=',num2str(out.montecarlo.p)])
    disp(' ')
    disp('-- Conjunction testing --')
    disp(['Indirect effect X->M->Y: P(ab=0|H0)=',num2str(out.conj.p)])
    disp(' ')
end
if display
    out.hf = figure('color',[1 1 1],'name','mediation analysis');
    
    % un-mediated model
    ha = subplot(4,1,1,'parent',out.hf,'nextplot','add');
    hc(1) = rectangle('parent',ha,'position',[-2 0 1 1],'curvature',[1 1]);
    ht(1) = text(-1.5,0.5,'X','parent',ha,'HorizontalAlignment','Center','FontSize',16);
    hc(2) = rectangle('parent',ha,'position',[4 0 1 1],'curvature',[1 1]);
    ht(2) = text(4.5,0.5,'Y','parent',ha,'HorizontalAlignment','Center','FontSize',16);
    % effet of X on Y in the absence of M
    if p1 > alpha
        col = 'r';
    else
        col = 'g';
    end
    hq(1) = quiver(-1,0.5,5.5,0,'parent',ha,'color',col);
    ht(4) = text(1.5,0.65,['R2=',num2str(round(1e3*a1.R2_a)/10),'% (p=',num2str(out.p(1,1)),')'],'parent',ha,'HorizontalAlignment','Center');
    axis(ha,'equal')
    axis(ha,'off')
    
    % mediated model
    ha = axes('parent',out.hf,'nextplot','add','position',[0.1300    0.1100    0.7750    0.5],'units','normalized');
    hc(1) = rectangle('parent',ha,'position',[-2 0 1 1],'curvature',[1 1]);
    ht(1) = text(-1.5,0.5,'X','parent',ha,'HorizontalAlignment','Center','FontSize',16);
    hc(2) = rectangle('parent',ha,'position',[4 0 1 1],'curvature',[1 1]);
    ht(2) = text(4.5,0.5,'Y','parent',ha,'HorizontalAlignment','Center','FontSize',16);
    hc(3) = rectangle('parent',ha,'position',[1 -2 1 1],'curvature',[1 1]);
    ht(3) = text(1.5,-1.5,'M','parent',ha,'HorizontalAlignment','Center','FontSize',16);
    % direct effet 
    if p3 > alpha
        col = 'r';
    else
        col = 'g';
    end
    hq(1) = quiver(-1,0.5,5.5,0,'parent',ha,'color',col);
    ht(4) = text(1.5,0.65,['R2=',num2str(round(1e3*de)/10),'% (p=',num2str(out.p(1,2)),')'],'parent',ha,'HorizontalAlignment','Center');
    % effect of X on M
    if p2 > alpha
        col = 'r';
    else
        col = 'g';
    end
    hq(2) = quiver(0.5*cos(pi/6)-1.5,0.5-0.5*sin(pi/6),2+0.5*cos(pi/6),sin(pi/6)-2,'parent',ha,'color',col);
    ht(5) = text(0,-0.5,['R2=',num2str(round(1e3*a2.R2_a)/10),'% (p=',num2str(out.p(2,1)),')'],'parent',ha,'HorizontalAlignment','Center');
    % effect of M on Y (in the presence of X)
    if p4 > alpha
        col = 'r';
    else
        col = 'g';
    end
    hq(3) = quiver(0.5*cos(pi/6)+1.5,0.5*sin(pi/6)-1.5,2+0.5*cos(pi/6),2-sin(pi/6),'parent',ha,'color',col);
    ht(5) = text(3,-0.5,['R2=',num2str(round(1e3*a4.R2_a)/10),'% (p=',num2str(out.p(2,2)),')'],'parent',ha,'HorizontalAlignment','Center');
    % indirect effect (Sobel test)
    ht(5) = text(1.5,-2.3,['Indirect effect [X->M->Y]: R2=',num2str(round(1e3*ie)/10),'% (Sobel test: p=',num2str(out.sobel.p),')'],'parent',ha,'HorizontalAlignment','Center');
    axis(ha,'equal')
    axis(ha,'off')
    
end


