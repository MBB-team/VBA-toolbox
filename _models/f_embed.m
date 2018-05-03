function [fx,J,dfdp] = f_embed(Xt,Theta,ut,in)
% evolution function for dynamical system's delay embedding
% function [fx,J,dfdp] = f_embed(Xt,Theta,ut,in)

fxd = zeros(in.dim.n,1);
Jd = zeros(in.dim.n_embed+in.dim.n,in.dim.n);
dfdpd = zeros(size(Theta,1),in.dim.n);
if sum(in.options.delays(:)) > 0
    for i=1:in.dim.n
        % First form delayed state vector and input
        iX = (1:in.dim.n) + (in.options.delays(i,:).*in.dim.n);
        dX = Xt(iX,:);
        % Evaluate evolution function at the delayed state vector:
        [opt,dim] = getOptions4EvalFun(in);
        [fx0,J0,dfdp0] = VBA_evalFun('f',dX,Theta,ut,opt,dim,0);
        fxd(i) = fx0(i);
        Jd(iX,i) = J0(:,i);
        dfdpd(:,i) = dfdp0(:,i);
    end
else
    % Evaluate evolution function at the non-delayed state vector:
    [opt,dim] = getOptions4EvalFun(in);
    [fxd,J0,dfdpd] = VBA_evalFun('f',Xt(1:in.dim.n,:),Theta,ut,opt,dim,0);
    Jd(1:in.dim.n,:) = J0;
end

% Construct full embedding flow and gradients:
Xe = Xt(1:in.dim.n_embed);
fx = [fxd;Xe];
J = zeros(in.dim.n_embed+in.dim.n);
J(:,1:in.dim.n) = Jd;
J(1:in.dim.n_embed,in.dim.n+1:in.dim.n+in.dim.n_embed) = eye(in.dim.n_embed);
dfdp = [dfdpd,zeros(in.dim.n_theta,in.dim.n_embed)];


function [opt,dim] = getOptions4EvalFun(in)
opt = in.options;
opt.f_fname = in.f_fname;
opt.f_nout = in.f_nout;
opt.checkGrads = 0;
dim.n = in.dim.n;
dim.n_theta = in.dim.n_theta;