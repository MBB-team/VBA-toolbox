function mov=prepare_DCMmovie(posterior,out)



%% loop
% get timeseries
[XS,dfdx,n_t,YS]=simulate_micro_u(out,posterior);

nodes_idx = out.options.inF.n5;
n_nodes = numel(nodes_idx);
inF = out.options.inF;


for t=1:n_t
   
mov(t).activity  = XS(nodes_idx,t);
connectivity = dfdx{t} / inF.deltat;
connectivity = connectivity(nodes_idx,nodes_idx);
mov(t).response  = YS(end-1:end,t);


mov(t).connectivity = connectivity -  diag(diag(connectivity));
mov(t).connectivity = mov(t).connectivity(:);
mov(t).pattern = get_pattern(out);
mov(t).time = (t-1)*out.options.inF.deltat;
mov(t).input = out.u(:,t)';

end

  % normalization
  max_act = max(abs([mov.activity]'))';
  max_con = max(abs([mov.connectivity]'))';
  
  for t=1:n_t
      mov(t).activity = mov(t).activity  ./ max(max_act);
      mov(t).connectivity = reshape(mov(t).connectivity(:) ,n_nodes,n_nodes)' ;%./ max(max_con)
      for i=1:n_nodes
        mov(t).connectivity(i,:) = mov(t).connectivity(i,:) .* mov(t).activity';
     end
  end
      

end

function [XS,dfdx,n_t,YS]=simulate_micro_u(out,posterior)
    n_t = out.dim.n_t*out.options.decim;
    Theta = posterior.muTheta;
    Phi = posterior.muPhi;
    US = out.u;
    inF = out.options.inF;
    inG = out.options.inG;
    
    X0 = posterior.muX0;
    XS = zeros(numel(X0),n_t);
    YS = zeros(8,n_t);
    XS(:,1) = feval(out.options.f_fname,X0,Theta,US(:,1),inF);
    
    for t=1:n_t
        [XS(:,t+1),dfdx{t},~] = feval(out.options.f_fname,XS(:,t),Theta,US(:,t),inF);
        YS(:,t) = feval(out.options.g_fname,XS(:,t),Phi,US(:,t),inG);
    end
    
end

function [pattern]=get_pattern(out)


    inF = out.options.inF;
    n_x = out.dim.n;%numel(inF.n5);
    n_u = out.dim.u;
    n_theta = out.dim.n_theta;
    base = ones(n_theta,1);
    idx=out.options.inF.n5;

    [~,to_keep,~]=f_DCMwHRFext(ones(n_x,1),base,ones(n_u,1),out.options.inF);
    pattern = to_keep(idx,idx)';

%     pattern = inF.A > 0;
%     
%     for i=1:n_u
%         pattern = pattern | (inF.B{i}>0);
%     end
%     for i=1:n_x
%         pattern = pattern | (inF.D{i}>0);
%     end
    
end