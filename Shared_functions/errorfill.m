function figurehandle= errorfill(x, y, varargin)

% Usage: figurehandle= errorfill(x, y, E1, E2, ..., LineSpec)
% 
% Use this function to draw a curve and its areas of confidence.
% Draws the line y(x), then draws the areas E1(x), E2(x), ...
% (Well actually it's done in the opposite order, to minimize overlap.)
% 
% x:        The x values for y, E1, E2, ... 
%           If x is a scalar, then (x, x+1, x+2, ...) is used.
%           If x is empty, then (1, 2, ...) is used.
%           Values must be finite and not NaN.
% 
% y:        Defines the curve, y(x).
%           Values must be finite and not NaN.
% 
% E:        The 'error' of y. Defines the area.
%           If E is a scalar, the shade will be between y(n)*(1+E) and y(n)*(1-E).
%           If E is of size (1,2) or (2,1), it will be between y(n)*(1+E(1)) and y(n)*(1-E(2)).
%           If E is a vector, the shade will be between y(n)+E(n) and y(n)-E(n).
%           If E is a double vector, the shade will be between y(n)+E(1,n) and y(n)-E(2,n).
%           Inf is marked by a jump up to a '^' symbol.
%           -Inf is marked by a jump down to a 'v' symbol.
%           NaN is marked by a jump up or down to a '*' symbol.
% 
% LineSpec: (Default: 'k') Line specification for y(x), e.g. 'b+'.
% 
% NOTE:     The E:s must be given in increasing size, or an area will be covered by a 
%           previous and larger area.
% 
% COMPATIBILITY:
%           errorbar(x, y, E, LineSpec)    corresponds to errorfill(x, y, E, LineSpec).
%           errorbar(x, y, L, U, LineSpec) corresponds to errorfill(x, y, [U; L], LineSpec).
% 
% EXAMPLE:  x= (1:0.25:10);    y= x.^2;
%           E3= [ones(1, length(x))*30; 0.7*y];  % 30 above and 0.3*y below
%           E3(1,7)= Inf;      E3(2,10)= Inf;       E3(2,25)= NaN;
%           errorfill(x, y, [0.1 0.02], 0.3, E3, 'b-+');
%           legend('E3', 'E2', 'E1', 'x^2','Inf','-Inf','NaN');
% 
% INSPIRATION: jackknife and confplot:
%           http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=10740
%           http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=2683
% 
% KEYWORDS: plot, errorplot, confidence interval, errors, patch, shaded plot, shading
% 
% TO DO:    Make it work for y-values of Inf, -Inf, and NaN.
% 
% Copyright (C) Peder Axensten (peder at axensten dot se), 2006.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version 1.0, submitted to Matlab Central File Exchange 2006-04-23:
% - First version in a presentable form.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	lSpec=		'k';												% Default values.
	N=			length(varargin);
	
	% Check arguments.
	if(ischar(varargin{end}))										% Get LineSpec.
		lSpec=		[lSpec varargin{end}];
		N=			N-1;
	end
	
	if(N < 1), error('Need more arguments!'); end					% Arguments?
	len=		length(y);											% Vector length.
	
	if(~isnumeric(y) || (len < 2) || (len ~= numel(y)))				% Check y.
		error('Argument 2 (y) must be a numeric vector.');
	end
	
	if(isempty(x))													% Check x.
		x=			1;
	elseif(~isnumeric(x))
		error('Argument 1 (x) must be numeric (scalar or vectror).');
	end
	
	if(length(x) ~= length(y))										% Check x & y.
		if(length(x) > 1)
			error('Argument 1 (x) and 2 (y) must have same size (or 1 must be a scalar).');
		else
			x=			(x:x+len-1);								% Make vector of scalar.
		end
	end
	
	% Get colors.
	c=			regexprep(lSpec, '([^rgbcmykw])', '');				% Calculate color.
	lSpec=		regexprep(lSpec, '([rgbcmykw])', '');				% (Remove [double] color.)
	cTable=		[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 1 1 1];
	color=		cTable('rgbcmykw' == c(end),:);				% Line color.
	colorshade=	reshape(1 + (1:N)'*(color-1)/N/2, 1,N,3);			% Calculate the shades.
	
	
	x=			transpose(repmat([x x(end:-1:1)], N, 1));			% x "there and back".
	hold on;														% Start the patches.
	
	% Get error values.
	for n= 1:N														% Check E.
		temp=	varargin{N-n+1};
		if(~isnumeric(temp) )
			error('Argument %d must be a non-empty numeric.', n+2);
		elseif(length(temp) == 1)									% A scalar.
			E(:,n)=		[y.*(1+temp)     y(end:-1:1).*(1-temp)];
		elseif(length(temp) == 2)									% A double scalar.
			E(:,n)=		[y*(1+temp(1))   y(end:-1:1)*(1-temp(2))];
		elseif(numel(temp) == len)									% A vector.
			E(:,n)=		[y+temp          y(end:-1:1)-temp(end:-1:1)];
		elseif(numel(temp) == 2*len)								% Double vector.
			E(:,n)=		[y+temp(1,:)     y(end:-1:1)-temp(2,end:-1:1)];
		else														% Something else...
			error('Strange size of argument %d.', n+2);
		end
	end
	
	% Handle NaN, Inf and -Inf.
	xInf=			x(E == Inf);									% Infinity markers.
	xnInf=			x(E == -Inf);									% -Infinity markers.
	xNaN=			x(isnan(E(1:len,:)));							% NaN markers.
	xnNaN=			x(isnan(E));									% NaN markers.
	finiteE=		E(isfinite(E));									% Find (ir)regularities.
	addition=		abs(max([finiteE(:)' y]))/10;					% Above/below max/min.
	maxval=			max([finiteE(:)' y])+addition;					% Represents Inf.
	minval=			min([finiteE(:)' y])-addition;					% Represents -Inf.
	E((E == Inf))=	maxval;
	E(E == -Inf)=	minval;
	E(isnan(E(1:len,:)))=	maxval;
	E(isnan(E))=	minval;
	
	% Draw.
	% Tried, but created bad eps (loop instead): patch(x, E, colorshade, 'LineStyle', 'none');
	for n= 1:N
		bb=patch(x(1:2*len), E(:,n), colorshade(1,n,:), 'LineStyle', 'none');	% Draw error fill.
        
        %%%%% Change la transparence
        
        set(bb,'FaceAlpha',0.3); %0.15
    end
	plot(x(1:len), y, lSpec, 'Color', color, 'linewidth', 2);		% Draw y(x).
	if(~isempty(xInf)),  plot(xInf,  maxval, '^k'); end				% Draw Inf points.
	if(~isempty(xnInf)), plot(xnInf, minval, 'vk'); end				% Draw -Inf points.
	if(~isempty(xNaN)),  plot(xNaN,  maxval, '*k'); end				% Draw high NaN points.
	if(~isempty(xnNaN)), plot(xnNaN, minval, '*k'); end				% Draw low NaN points.
	figurehandle=	get(0, 'CurrentFigure');						% Return figure handle.
	hold off;														% We are done.
	
	% Notify.
	if(numel(xInf)+numel(xnInf)+numel(xNaN)+numel(xnNaN) > 0)		% Irregular values?
		warning('At least some arguments contain irregular values (see ''help errorfill'').');
	end
end
