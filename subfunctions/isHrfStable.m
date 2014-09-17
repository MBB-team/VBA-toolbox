function [flag] = isHrfStable(X,kas,kaf,phi)

try
    phi;
catch
    phi = -0.5785; 
end

% get vasodilatory signal and blood inflow
X1 = X(1);
X2 = X(2);

% if X(1) > 0
%     flag = 1;
%     return
% end

% first compute separatrix
x2 = [-5:1e-3:5]';
n = length(x2);
x1 = phi.*exp(kaf/2-abs(kas-kaf/2)).*exp(x2);
dx1 = diff(x1);
dx2 = diff(x2);

% get closest point on separatrix
d2 = sum([(x1-X1).^2,(x2-X2).^2],2);
[md,ind] = min(d2);

% see whether the point is inside the domain
vec = [x2(ind)-X2,x1(ind)-X1];
tan = [-dx1(ind),dx2(ind)];
nv = sqrt(sum(vec.^2));
nt = sqrt(sum(tan.^2));
angle = acos(sum(vec.*tan)./(nv*nt));
flag = -sign(cos(angle));

hf = figure;
ha = axes(...
    'parent',hf,...
    'nextplot','add');
plot(ha,x2,x1,'.')
plot(ha,x2(ind),x1(ind),'ro')
quiver(x2(1:end-1),...
    x1(1:end-1),...
    -diff(x1),...
    diff(x2),...
    2)
quiver(X2,...
    X1,...
    vec(1),...
    vec(2),...
	'color',[1 0 0])
grid(ha,'on')
axis(ha,'tight')

plot(ha,X2,X1,'r+')
plot(ha,X2,X1,'ro')


