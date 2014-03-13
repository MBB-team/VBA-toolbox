function connect = grapher_connectivityPattern(out,theta)


n_x=out.dim.n;
n_theta=out.dim.n_theta;
n_u=size(out.u,1);
idx=out.options.inF.n5;

if nargin<2
    theta = zeros(n_theta,1);
end

base = ones(n_theta,1);
[~,to_keep,~]=f_DCMwHRFext(ones(n_x,1),base,ones(n_u,1),out.options.inF);
to_keep = to_keep(idx,idx);

[~,connect,~]= f_DCMwHRFext(ones(n_x,1),theta(:),ones(n_u,1),out.options.inF);
connect = connect(idx,idx);
connect(to_keep==0) = NaN;
connect(1:numel(idx)+1:end)=NaN;
connect = connect /  out.options.inF.deltat;


end