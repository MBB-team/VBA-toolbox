function [r,t,p]=spear(x,y)
%Syntax: [r,t,p]=spear(x,y)
%__________________________
%
% Spearman's rank correalation coefficient.
%
% r is the Spearman's rank correlation coefficient.
% t is the t-ratio of r.
% p is the corresponding p-value.
% x is the first data series (column).
% y is the second data series, a matrix which may contain one or multiple
%     columns.
%
%
% Reference:
% Press W. H., Teukolsky S. A., Vetterling W. T., Flannery B. P.(1996):
% Numerical Recipes in C, Cambridge University Press. Page 641.
%
%
% Example:
% x = [1 2 3 3 3]';
% y = [1 2 2 4 3; rand(1,5)]';
% [r,t,p] = spear(x,y)
%
%
% Products Required:
% Statistics Toolbox
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% 45110- Dourouti
% Ioannina
% Greece
% 
% University e-mail: me00743@cc.uoi.gr
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
% 
% 3 Feb 2002.


% x and y must have equal number of rows
if size(x,1)~=size(y,1)
    error('x and y must have equal number of rows.');
end


% Find the data length
N = length(x);

% Get the ranks of x
R = crank(x)';

for i=1:size(y,2)
    
    % Get the ranks of y
    S = crank(y(:,i))';
    
    % Calculate the correlation coefficient
    r(i) = 1-6*sum((R-S).^2)/N/(N^2-1);
    
end

% Calculate the t statistic
if r == 1 | r == -1
    t = r*inf;
else
    t=r.*sqrt((N-2)./(1-r.^2));
end

% Calculate the p-values
p=2*(1-VBA_spm_Tcdf(abs(t),N-2));





function r=crank(x)
%Syntax: r=crank(x)
%__________________
%
% Assigns ranks on a data series x. 
%
% r is the vector of the ranks
% x is the data series. It must be sorted.
%
%
% Reference:
% Press W. H., Teukolsky S. A., Vetterling W. T., Flannery B. P.(1996):
% Numerical Recipes in C, Cambridge University Press. Page 642.
%
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% 45110- Dourouti
% Ioannina
% Greece
% 
% University e-mail: me00743@cc.uoi.gr
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
% 
% 3 Feb 2002.


u = unique(x);
[xs,z1] = sort(x);
[z1,z2] = sort(z1);
r = (1:length(x))';
r=r(z2);

for i=1:length(u)
    
    s=find(u(i)==x);
    
    r(s,1) = mean(r(s));
    
end




