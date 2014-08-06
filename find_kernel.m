function [kernel,y] = find_kernel(posterior,out,nt)

if nargin < 3
    nt = 8;
end
u=out.u;
nu = size(u,1);

out.options.verbose = 0;
out.options.DisplayWin = 0;
% 
% if out.options.microU == 1
%     [~,gx,microTime,sampleInd] = VBA_microTime(posterior,u,out) ;
%     y = gx(out.options.sources(2).out,2:end);
% else
%     [yp,~,~,~,er] = simulateNLSS(...
%         out.options.dim.n_t,...
%         out.options.f_fname,...
%         out.options.g_fname,...
%         posterior.muTheta,...
%         posterior.muPhi,...
%         out.u,...
%         Inf,Inf,out.options,...
%         posterior.muX0);
%     resp_idx = [out.options.sources(2:end).out];
%     y = yp(resp_idx,:)-er(resp_idx,:);
% end
y=out.suffStat.gx(out.options.sources(2).out,:);

% kernel=[];
% return;
%%%%%%%



% if out.options.microU == 1
%     for iu=1:nu
%     u_macro(iu,:) = mean(reshape(u(iu,:),out.options.decim,round(size(u,2)/out.options.decim)),1);
%     end
%     u=u_macro;
% end

%     for iu=1:nu
%     u_micro(iu,:) = vec(repmat(u(iu,:),out.options.decim,1));
%     end
%     u=u_micro;

% shorten
n_t=size(y,2);
% u=u([1 2 5 6],1:n_t);
y=y(:,1:n_t);
% nu=4;

inG.dim.n_t = nt;
inG.deltat = out.options.inF.deltat*out.options.decim;
inG.dim.nu = nu;
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = 3*inG.dim.nu +1;

opt.inG = inG;
opt.priors.muPhi = [repmat([0;.1;1],inG.dim.nu,1);0]; 
opt.priors.SigmaPhi = 10*kron(eye(inG.dim.nu),diag([10,5,10]));
opt.priors.SigmaPhi(end+1, end+1) = 1;
opt.verbose = 1;
opt.DisplayWin = 1;
opt.checkGrads = 0;
opt.TolFun = 1e-8;
opt.GnTolFun = 1e-8;
opt.MaxIter = 100;
opt.GnMaxIter = 100;

%normalize u for repetitions
% for i=1:nu
%    cpt = [0];
%    last = u(i,1);
%    for t=2:size(u,2)
%        if u(i,t) ~= last
%            if last == 0
%                cpt(end+1) = 1;
%            end
%        else
%            if last ~= 0
%                cpt(end) = cpt(end)+1;
%            end
%        end 
%        last = u(i,t);
%    end
%    u(i,:) = u(i,:)/max(cpt);
% end
u(2,:) = min(u(2,:),u(1,:));

kernel = struct();
for iObs = 1:size(y,1)
    g_fname = @g_convSig_approx;
    

[pk,ok] = VBA_NLStateSpaceModel(vec(y(iObs,:)) ,vec(u'),[],g_fname,dim,opt);

    kernel(iObs).params = reshape(pk.muPhi(1:end-1),3,nu);
    for iu = 1:nu
      % first moment 
%       [kernel(iObs).timeseries(iu,:),~,kernel(iObs).landmarks(iu,:),dadp] = kernel_sinexp(kernel(iObs).params(1,iu),kernel(iObs).params(2,iu),kernel(iObs).params(3,iu),(0:nt-1)*inG.deltat) ;
      [kernel(iObs).timeseries(iu,:),~,kernel(iObs).landmarks(iu,:),dadp] = kernel_sinexp(kernel(iObs).params(1,iu),kernel(iObs).params(2,iu),kernel(iObs).params(3,iu),0:.01:(nt-1)*inG.deltat) ;
      % second moment
      idx = (iu-1)*3 + (1:3) ;
      kernel(iObs).sigma.landmarks(iu) = dadp'*pk.SigmaPhi(idx,idx)*dadp;
    end
    kernel(iObs).sigma.params = reshape(diag(pk.SigmaPhi(1:end-1,1:end-1)),3,nu);
    
   
close all; plot(ok.suffStat.gx,'r'); hold on; plot(ok.y); hold off  
pause
    
end

