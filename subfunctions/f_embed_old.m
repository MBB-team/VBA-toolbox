function [fx,J,dfdp] = f_embed(Xt,Theta,ut,in)
% evolution function for dynamical system's dimension embedding
% function [fx,J,dfdp] = f_embed(Xt,Theta,ut,in)

% First evaluate standard evolution function
P = Theta(1:in.dim.n_theta);
Xt0 = Xt(1:in.dim.n);
[opt,dim] = getOptions4EvalFun(in);
[fx0,J0,dfdp0] = VBA_evalFun('f',Xt0,P,ut,opt,dim);

% Then evaluate additional hidden states flow and derivatives
Xe = Xt(in.dim.n+1:in.dim.n+in.dim.n_embed);
In = eye(in.dim.n);
[SP,dsdp] = sigm(Theta(in.dim.n_theta+1:end),...
    struct('G0',2,'S0',-1,'beta',1,'INV',0));
Je = reshape(SP,in.dim.n,in.dim.n_embed);
fxe = in.dt.*Je*Xe;
dfdpe = diag(dsdp)*kron(Xe,In);
Zeros = zeros(in.dim.n,in.dim.n_embed);
fx      = [ fx0 + fxe
            Xe         ];
J       = [ J0          Zeros
            in.dt*Je'   eye(in.dim.n_embed) ];
dfdp    = [ dfdp0       zeros(size(dfdp0,1),in.dim.n_embed)
            in.dt*dfdpe zeros(size(dfdpe,1),in.dim.n_embed) ];
       
        
function [opt,dim] = getOptions4EvalFun(in)
opt = in.options;
opt.f_fname = in.f_fname;
opt.f_nout = in.f_nout;
opt.decim = 1;
dim.n = in.dim.n;
dim.n_theta = in.dim.n_theta;