function gplotdc(A,xy,varargin)
%GPLOTDC Plot a Directed Graph
% GPLOTDC(A,XY) Plots the Directed Graph represented by adjacency
%   matrix A and points xy using the default style described below
% GPLOTDC(A,XY,PARAM1,VAL1,PARAM2,VAL2,...) Plots the Directed Graph
%   using valid parameter name/value pairs
%
%   Inputs:
%       A     - NxN adjacency matrix, where A(I,J) is nonzero (=1)
%               if and only if there is an edge between points I and J
%       xy    - Nx2 matrix of x/y coordinates
%       ...   - Parameter name/value pairs that are consistent with
%               valid PLOT parameters can also be specified
%
%   Default Plot Style Details:
%       1. Undirected (2-way) edges are plotted as straight solid lines
%       2. Directed (1-way) edges are plotted as curved dotted lines with
%           the curvature bending counterclockwise moving away from a point
%       3. Any vertex that is connected to itself is plotted with a
%           circle around it
%
%   Example:
%       % plot a directed graph using default line styles
%       n = 9; t = 2*pi/n*(0:n-1);
%       A = round(rand(n));
%       xy = [cos(t); sin(t)]';
%       gplotdc(A,xy);
%       for k = 1:n
%           text(xy(k,1),xy(k,2),['  ' num2str(k)],'Color','k', ...
%               'FontSize',12,'FontWeight','b')
%       end
%
%   Example:
%       % plot a directed graph using plot name/value parameter pairs
%       n = 9; t = 2*pi/n*(0:n-1);
%       A = round(rand(n));
%       xy = [cos(t); sin(t)]';
%       gplotdc(A,xy,'LineWidth',2,'MarkerSize',8);
%
%   Example:
%       % plot a directed graph using a different color and linestyle
%       n = 9; t = 2*pi/n*(0:n-1);
%       A = round(rand(n));
%       xy = [cos(t); sin(t)]';
%       gplotdc(A,xy,'Color',[0.67 0 1],'LineStyle',':');
%
% See also: gplot, plot
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Release Date: 4/12/08

% Process Inputs
if nargin < 2
    error('Not enough input arguments.');
end
[nr,nc] = size(A);
[n,dim] = size(xy);
if (~n) || (nr ~= n) || (nc ~= n) || (dim < 2)
    eval(['help ' mfilename]);
    error('Invalid input. See help notes above.');
end
params = struct();
for var = 1:2:length(varargin)-1
    params.(varargin{var}) = varargin{var+1};
end

% Parse the Adjacency Matrix
A = double(logical(A));
iA = diag(diag(A));         % self-connecting edges
dA = A.*(1-A');             % directed edges (1-way)
uA = A-dA-iA;               % undirected edges (2-way)

% Make NaN-Separated XY Vectors
[ix,iy] = makeXY(iA,xy);
[dx,dy] = makeXY(dA,xy);
[ux,uy] = makeXY(tril(uA,0),xy);

% Add Curvature to Directed Edges
dX = dx;
dY = dy;
pct = 0.04;
for k = 1:4
    [dX,dY,pct] = makeCurved(dX,dY,pct);
end

% Plot the Graph
plot(ux,uy,'b-',params)
hold on
plot(dX,dY,'r--',params)
plot(ix,iy,'ko',params)
plot(xy(:,1),xy(:,2),'k.')
hold off

    function [x,y] = makeXY(A,xy)
        if any(A(:))
            [J,I] = find(A');
            m = length(I);
            xmat = [xy(I,1) xy(J,1) NaN(m,1)]';
            ymat = [xy(I,2) xy(J,2) NaN(m,1)]';
            x = xmat(:);
            y = ymat(:);
        else
            x = NaN;
            y = NaN;
        end
    end

    function [X,Y,PCT] = makeCurved(x,y,pct)
        N = length(x);
        if N < 2
            X = x;
            Y = y;
        else
            M = 2*N-1;
            X = zeros(1,M);
            Y = zeros(1,M);
            X(1:2:M) = x;
            Y(1:2:M) = y;
            X(2:2:M-1) = 0.5*(x(1:N-1)+x(2:N))+pct*diff(y);
            Y(2:2:M-1) = 0.5*(y(1:N-1)+y(2:N))-pct*diff(x);
        end
        PCT = 0.5*pct;
    end
end
