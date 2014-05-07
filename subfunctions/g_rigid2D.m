function [gx,dgdx,dgdp] = g_rigid2D(x,P,u,in)
% apply rigid body 2D spatial transform (rotation and translation)
% function [gx,dgdx,dgdp] = g_rigid2D(x,P,u,in)
% IN:
%   - x: [useless]
%   - P: the 3 rigid body transform parameters
%   - u: [useless]
%   - in: contains the original image (in.Y)
% OUT:
%   - gx: vectorized transformed image
%   - dgdx: []
%   - dgdp: gradient of the image wrt to transform parameters

gx = applyRT(P,in);

% the numerical derivation is done here because of interpolation errors:
dgdx = [];
dgdp = zeros(size(gx,1),3);
I = eye(3);
for i=1:3
    dgdp(:,i) = applyRT(P+I(:,i),in) - gx;
end
dgdp = dgdp';


function gx = applyRT(P,in)

tx = P(1);
ty = P(2);
th = P(3);
T = [ 1 0 tx
      0 1 ty
      0 0 1  ];
R = [ cos(th) -sin(th) 0
      sin(th) cos(th)  0
      0       0        1  ];
[nx,ny] = size(in.Y);
[xg,yg] = meshgrid(-nx/2+1:nx/2,-ny/2+1:ny/2);
Txy = T*[xg(:)';yg(:)';ones(1,nx*ny)];
RTxy = R*[Txy(1,:);Txy(2,:);ones(1,nx*ny)];
RTx = reshape(RTxy(1,:)',nx,ny);
RTy = reshape(RTxy(2,:)',nx,ny);
gx = griddata(RTx,RTy,in.Y,xg,yg,'cubic');
gx = gx(:);
gx(isnan(gx)) = 0;

