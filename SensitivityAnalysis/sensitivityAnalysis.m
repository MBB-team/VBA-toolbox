function [results]=sensitivityAnalysis(posterior,out,obs,contrast)




theta = posterior.muTheta;
phi = posterior.muPhi;

options=out.options;
dim=options.dim;

u=out.u(:,1:dim.n_t);

if nargin <3
    contrast=out.u;
end
contrast=contrast(:,1:dim.n_t);

results.contrast=contrast;

%== get parameters of interest
inF= out.options.inF;
thetas=[inF.indA inF.indB{:} inF.indC inF.indD{:}];

thetas_lbl={};
for i=1:length(inF.indA)
    thetas_lbl{end+1} = sprintf('A_%d',i);
end
for k=1:length(inF.indB)
    for i=1:length(inF.indB{k})
    thetas_lbl{end+1} = sprintf('B^%d_%d',k,i);
    end
end
for i=1:length(inF.indC)
    thetas_lbl{end+1} = sprintf('C_%d',i);
end
for k=1:length(inF.indD)
    for i=1:length(inF.indD{k})
    thetas_lbl{end+1} = sprintf('D^%d_%d',k,i);
    end
end

results.labels=thetas_lbl;


%== simulate trajectory

%- initialization
x = zeros(dim.n,dim.n_t);
y = zeros(dim.p,dim.n_t);
dgdp = zeros(dim.n_theta,dim.p,dim.n_t);


x(:,1) = posterior.muX0;
y(:,1) = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim,1);

%- loop over time
for t = 2:dim.n_t
        
    % Evaluate evolution function at past hidden state
    if dim.n > 0 
        [x(:,t),~,dfdp] = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim,t);
    end
    
    % Evaluate observation function at current hidden state
    [y(:,t),dgdx,~] = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim,t);
    
    % compute gradient
    dgdp(:,:,t)= dfdp*dgdx;


end

%== summarize 
%- get coordinate system in input space
perturb=SplitGrad(dgdp(thetas,:,:),contrast);
predict=SplitGrad(reshape(y,[1,size(y)]),contrast);
effect=SplitGrad(reshape(out.y,[1,size(out.y)]),contrast);
SigmaEffect=SplitGrad(reshape(out.suffStat.vy,[1,size(out.y)]),contrast);

results.effect=effect;
results.SigmaEffect=SigmaEffect;
results.predict=predict;
results.perturb=perturb;

%== decompose effects
% for i=1:size(perturb,1)
%     base = cell2mat(cellfun(@(x) x(:)/norm(x(:)),{perturb{i,:}},'UniformOutput',false)) ;%/norm(x(:))
%     
% %     base = [base eye(size(base,1))] %/size(base,1)
%     base(:,end+1)=1/sqrt(size(base,1))
%     w1 =  base\(predict{i}(:));
% 
%     res{i} = predict{i}(:)-base*w1(:) ;
%     w{i} = w1 / norm(w1(1:end-1)); %(1:end-size(base,1))
%     
%     
% end

% get sources types for offset
% stypes=[];
% for i=1:numel(out.options.sources)
%   stypes(out.options.sources(i).out) =  out.options.sources(i).type;
% end

for i=1:size(perturb,1)
    
    for j=1:length(thetas)
       base = perturb{i,j}(:) ;
       base(:,end+1)=1;
       w1 =  base\(predict{i}(:));
       res = predict{i}(:)-base*w1(:) ;
%        w2(j) = 1/(norm(res).^2*norm(w1(1))^2); % lat>all- for>all(+for)
%        w2(j) = 1/(norm(res).^2*norm(w1(1))); % lat>lat for>all(+lat)
%         w2(j) = 1/(norm(res).^2*norm(w1)^2); % lat>all- for>all(+for)
       w2(j) = 1/(norm(res).^2*norm(w1)); % lat>lat for>all(+lat)
       
       
        
    end
    % check for zero scored
     goodw = isfinite(w2);
     w2 = w2(:)/norm(w2(goodw));
     w2(~goodw) = 0;
     w{i} = w2;
end


results.w = w;
% results.residuals= cellfun(@reshape,res(:),cellfun(@size,effect,'UniformOutput',false), 'UniformOutput',false);

%==
sensitivityDisplay(results);
end

function r=SplitGrad(grad,u)
for nu=1:size(u,1)
   base_u{nu} = unique(u(nu,:));  
   
   base_u{nu} = setdiff(base_u{nu},[0]);
end
if length(base_u)>1
    categs=product(base_u);
else
    categs = base_u{1};
end
%normalization

for y_i=1:size(grad,2)
    for theta_i=1:size(grad,1)
        dim_r = cellfun(@length,base_u);
%         if length(dim_r)==1;
%             dim_r = [1 dim_r];
%         end
    r{y_i,theta_i} = zeros([dim_r 1]);
    for i=1:size(categs,2)
        level=categs(:,i);
        idx = find(sum(repmat(level,1,size(u,2))==u,1)==size(categs,1));
        if ~isempty(idx)
            r{y_i,theta_i}(i)=mean(grad(theta_i,y_i,idx),3);
        else
            r{y_i,theta_i}(i)=NaN;
        end
    end
    end
end
end

% function y=getTraj(options,dim,x0,theta,phi,u)
%     x = zeros(dim.n,dim.n_t);
%     y = zeros(dim.p,dim.n_t);
% 
%     x(:,1) = x0;
%     y(:,1) = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim,1);
% 
% %- loop over time
% for t = 2:dim.n_t
%         
%     % Evaluate evolution function at past hidden state
%     if dim.n > 0 
%         x(:,t) = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim,t);
%     end
%     
%     % Evaluate observation function at current hidden state
%     y(:,t) = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim,t);
% end
% 
% end
% 
% function dgdp=getPerturb(options,dim,x0,theta,phi,u)
%     dgdp = zeros(dim.n_theta,dim.p,dim.n_t);
%     dthetha = 1e-6;
%     
%     for i=1:length(theta)
%         theta_p = zeros(length(theta),1);
%         theta_p(i) = dthetha;
%         yp = getTraj(options,dim,x0,theta+theta_p,phi,u) ;
%         ym = getTraj(options,dim,x0,theta-theta_p,phi,u) ;
%         dgdp_t = (yp-ym)/(2*dthetha);
%         dgdp(i,:,:)=reshape(dgdp_t,[1,dim.p,dim.n_t]);     
%     end
% 
%     
% 
% end



function P = product(sets)


% check for empty inputs
q = ~cellfun('isempty',sets) ;


ni = sum(q) ;
ii = ni:-1:1 ;

if ni==0,
    P = [] ;
else
    args = sets(q) ;
    if ni==1,
        P = args{1}(:) ;
    else
        % flip using ii if last column is changing fastest
        [P{ii}] = ndgrid(args{ii}) ;
        % concatenate
        P = reshape(cat(ni+1,P{:}),[],ni) ;
    end
end

P=P';
    
end

    


