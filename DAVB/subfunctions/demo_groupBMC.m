% this demo checks the behaviour group BMC

clear variables
close all


K = 8; % # models
n = 32; % # subjects

N = 16;
ep1 = zeros(K,N);
ep2 = zeros(K,N);

options.TolFun = 1e-4;
options.DisplayWin = 1;
options2 = options;
% for i=1:K
%     options.families{i} = i;
% end
options.families{1} = 1:floor(K/2);
options.families{2} = floor(K/2)+1:K;%floor(K/2);
% options.families{3} = floor(K/2)+1:floor(3*K/4);
% options.families{4} = floor(3*K/4)+1:K;

% hf = figure;
% ha(1) = subplot(2,2,1,'parent',hf,'nextplot','add');
% ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add');
% ha(3) = subplot(2,2,3,'parent',hf,'nextplot','add');
% 
% plot(ha(1),[0 1],[0 1],'r')
% title(ha(1),'model frequencies')
% plot(ha(2),[0 1],[0 1],'r')
% title(ha(2),'exceedance probabilities')
% plot(ha(3),[0 1],[0 1],'r')
% title(ha(3),'model attributions')

nf = length(options.families);
C = zeros(K,nf);
for i=1:nf
    indf = options.families{i};
    C(indf,i) = 1;
end


for ii=1:N
    
    ii
    
    bias = [10*ones(floor(K/2),1);zeros(ceil(K/2),1)];
    for i=1:n
        L(:,i) = 1*randn(K,1) + bias;
    end

    [posterior,out] = VBA_groupBMC(L,options);
    
    
    Lf = log(getFamily(L,options.families));
    
    [p2,o2] = VBA_groupBMC(Lf,options2);
    
%     pause
    % VBA_displayGroupBMC(posterior,out)
    
%     alpha0 = out.options.priors.a';
%     [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs(L',alpha0,1e4);
%     
%     ep1(:,ii) = out.ep;
%     ep2(:,ii) = xp;
%     
%     plot(ha(1),exp_r(:),out.Ef(:),'k.');
%     plot(ha(2),xp(:),out.ep(:),'k.');
%     g = g_post';
%     plot(ha(3),g(:),posterior.r(:),'k.');
%     
%     try
%         plot(out.options.handles.ha(3),exp_r,'go')
%         plot(out.options.handles.ha(4),xp,'go')
%     end
%     
%     drawnow
%     pause(1)
%     try,close(out.handles.hf);end
    
end





