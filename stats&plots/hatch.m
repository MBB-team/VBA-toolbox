function  hatch(obj,angle,color,style,step,width)
% HATCH  Hatches a two-dimensional domain.
%     Somewhat similar to the FILL command but fills the closed
%     contour with hatched lines instead of uniform color.
%     The boundary line(s) must be created prior to this command
%     and its handle(s) specified by OBJ argument.
%   HATCH(OBJ) Hatchs the domain bounded by the Xdata and Ydata
%     properties of the OBJ (can be line or patch type).
%   HATCH by itself takes as OBJ the last object created in the 
%     current axes.
%   HATCH(OBJ,ANGLE,COLOR,STYLE,STEP,WIDTH) Specifies additional
%     parameters:  slope of the hatches (in degrees),
%     their color ([red green blue], or 'r','g','b','w','y','c','m')
%     the linestyle ('-','--','-.',':') and also the steps
%     (distances between hatches) and the linewidth (thickness)
%     of the hatch lines (the last two in points).
%     All arguments are optional but must be in the given order.
%     They also can be grouped in two vectors for convenience:
%   HATCH(OBJ,[ANGLE STEP WIDTH],[COLOR STYLE]), where the last 
%     argument is a string: 'w--', '-.y' are both legal notations.
%   LH = HATCH(OBJ) also returns the handle of a hatch line.
%     OBJ can be a vector of several handles, in which case the
%     composite boundary is hatched. If one contour lies within
%     another, then the first one will appear as a "hole" in the 
%     outer contour.
%
%     Examples:
%   HATCH(L,30,[1 0 0],'--',8,2) or
%   HATCH(L,[30 8 2],'r--')    Hatchs a domain
%     bounded by the contour L with red dashed lines of 2-point
%     thickness and 8-point steps, inclined at 30 degrees
%     to the x-axis.
%   HATCH([L1 L2],'CROSS')  Hatches the domain with composite
%     boundary [L1 L2] with hatch parameters specified by the
%     "macro" 'CROSS'  (in this case two crossed hatches).
%
%	Type HATCH('demo') to see all the macro effects available
%
%     See also FILL, LINE, PATCH

%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  April 27 1994
%
%  Modified for Matlab 5+, Iram Weinstein weinsteini@saic.com
%  May 10, 2001



% Defaults ............
angledflt = 45;        % Angle in degrees
colordflt = [1 1 1];   % Color
styledflt = '-';       % Linestyle
widthdflt = 1;         % Thickness of the lines
stepdflt = 10;         % Distance between hatches
widthdflt = 1;         % Thickness of the lines

%macros .....................
macros={ 'hor'  '0.1	,''w''';
	'ver'  '90.1,''w''';
	'thick'     '[45 10 2],''w''';     % Thicker hatch lines
	'thicky'    '[45 10 2],''y''';
	'Thick'     '[45 10 2],''w''';
	'dense'     '[45 5]';            % Denser lines
	'Dense'     '[45 3]';
	'fill'      '[45 1 1]';        % Fills almost uniformly
	'filly'     '[45 1 1],''y''';
	'rare'      '[45 15]';
	'Rare'      '[45 20]';
	'rarethick' '[45 15 2]';
	'cross'     '45 135';            % Two crossed hatches
	'plus'      '1 90.1'};
argcount=nargin;
% Handle input ::::::::::::::::::::::::::::::::::::::::::::::::::
% Check for macros ................
ismacro = 0;
if argcount==2
	if isstr(angle)
		macro = angle; ismacro = 1; angle = angledflt;
	end
elseif argcount == 1
	if isstr(obj), 
		if strcmp(obj,'demo'),
			demo;
			return
		else
			macro = obj; ismacro = 1; argcount = 0; 
		end
	end
end

if argcount==0   % Find the object on figure
	ch = get(gca,'child');
	if ~isempty(ch), obj = ch(1);
	else
		disp([10 '  Error: no object to hatch is found ' 10])
		return
	end
end

%  Hatch with macros ..........................................
if ismacro
	match=strcmp(macros(:,1),macro);
	if ~any(match),
		disp([10 '  Macro is not found' 10]), return
	end
	call=macros(match,:);
	call=call{2};
	isbr = cumsum((call=='[')-(call==']'));
	lc = length(call);
	n0 = [0 find(call==' '&isbr==0) lc];
	for jf = 2:length(n0)
		if nargout == 0, out = '';
		elseif nargout==1, out = ['lh(' num2str(jf-1) ')='];
		end
		str = [out 'hatch(obj,' call(n0(jf-1)+1:n0(jf)) ');'];
		eval(str)
	end
	return
end

if argcount<6, width = widthdflt; end
if argcount<5, step = stepdflt;   end
if argcount<4, style = styledflt; end
if argcount<3, color = colordflt; end
if argcount<2, angle = angledflt; end

% Check for step and width in one vector
if length(angle)>1
	step = angle(2);
	if length(angle)>2, width = angle(3); end
	angle = angle(1);
end
% Check for color and style in one string
if isstr(color)
	A = color(ones(8,1),:)==setstr(ones(length(color),1)*'wyrgbcmk')';
	n0 = find(any(A));
	str = color(any(A)==0);
	if ~isempty(n0), color = color(n0); end
	if ~isempty(str)
		iss = strcmp(str,'-')|strcmp(str,'--')|strcmp(str,'-.');
		if iss|strcmp(str,':'), style = str; end
	end
end

% Check for the object to be line or patch
typ = get(obj(1),'type');
if ~(strcmp(typ,'line')|strcmp(typ,'patch'))
	disp([10 '  Error: object must be either line or patch ' 10])
	return
end

angle = angle*pi/180;             % Degrees to radians
% Combine all objects into a single contour
x = []; y = [];
for jo = 1:length(obj)% Get x,y
	xx=get(obj(jo),'xdata');
	x = [x xx(:).'];
	yy=get(obj(jo),'ydata');
	y = [y yy(:).'];
	if jo == 1,
		yi = find(~isnan(x)&~isnan(y));
		if ~isempty(yi), n0 = yi(1); x0 = x(n0); y0 = y(n0);
		else, x0 = 0; y0 = 0;
		end
	end
	x = [x x0]; y = [y y0];         % Close loop
end
yi = find(~isnan(x)&~isnan(y));
x = x(yi); y = y(yi);                       % Remove NaN's
ll = length(x);

% Transform the coordinates .............................
oldu = get(gca,'units');
set(gca,'units','points')
sza = get(gca,'pos'); sza = sza(3:4);
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
islx = strcmp(get(gca,'xscale'),'log');
isly = strcmp(get(gca,'yscale'),'log');
if islx     % If log scale in x
	xlim = log10(xlim);
	x = log10(x);
end
if isly     % If log scale in y
	ylim = log10(ylim);
	y = log10(y);
end
xsc = sza(1)/(xlim(2)-xlim(1)+eps);
ysc = sza(2)/(ylim(2)-ylim(1)+eps);

ca = cos(angle); sa = sin(angle);
x0 = mean(x); y0 = mean(y);  % Central point
x = (x-x0)*xsc; y = (y-y0)*ysc;
yi = x*ca+y*sa;              % Rotation
y = -x*sa+y*ca;
x = yi;
y = y/step;    % Make steps equal to one

% Compute the coordinates of the hatch line ...............
yi = ceil(y);
ll = length(y);
yd = [yi(2:ll)-yi(1:ll-1) 0];
dm = max(abs(yd));
fnd = find(yd);
lfnd = length(fnd);
A = sign(yd(fnd));
edm = ones(dm,1);
A = A(edm,:);
if size(A,1)>1, A = cumsum(A); end
fnd1 = find(abs(A)<=abs(yd(edm,fnd)));
A = A+yi(edm,fnd)-(A>0);
xy = (x(fnd+1)-x(fnd))./(y(fnd+1)-y(fnd));
xi = x(edm,fnd)+(A-y(edm,fnd)).*xy(edm,:);
yi = A(fnd1);
xi = xi(fnd1);

% Sorting points of the hatch line ........................
li = length(xi); 
xi0 = min(xi); xi1 = max(xi);
yi0 = min(yi); yi1 = max(yi);
ci = yi*(xi1-xi0)+xi;
[ci,num] = sort(ci);
xi = xi(num); yi = yi(num);
if floor(li/2)~=li/2
	xi = [xi xi(li)];
	yi = [yi yi(li)];
end

% Organize to pairs and separate by  NaN's ................
li = length(xi);
xi = reshape(xi,2,li/2);
yi = reshape(yi,2,li/2);
xi = [xi; ones(1,li/2)*nan];
yi = [yi; ones(1,li/2)*nan];
xi = xi(:)'; yi = yi(:)';

% Transform to the original coordinates ...................
yi = yi*step;
xy = xi*ca-yi*sa;
yi = xi*sa+yi*ca;
xi = xy/xsc+x0;
yi = yi/ysc+y0;
if islx, xi = 10.^xi; end
if isly, yi = 10.^yi; end

% Now create a line to hatch ..............................
ax=axis;
lh = line('xdata',xi,'ydata',yi);
set(lh,'linewidth',width);
set(lh,'color',color)
set(lh,'linestyle',style)
set(gca,'units',oldu)   % Set axes units back
axis(ax);
%------------------------
function demo
figure
subplot 211
h=bar(ones(2,7))
axis([0.6 1.4 0 1.3])
set(gca,'color',[1 1 1]*0.9,'XTickLabel','','YTickLabel','')
hatch(h(1),'hor')
hatch(h(2),'ver')
hatch(h(3),'thick')
hatch(h(4),'thicky')
hatch(h(5),'Thick')
hatch(h(6),'dense')
hatch(h(7),'Dense')

x=0.65+[0 1 2 3 4 5 6]*.11;
y=1.1*ones(size(x));
string=char({'hor', 'ver', 'thick','thicky','Thick','dense','Dense'});
text(x,y,string);

subplot 212
h=bar(ones(2,7))
axis([0.6 1.4 0 1.3])
set(gca,'color',[1 1 1]*0.9,'XTickLabel','','YTickLabel','')
hatch(h(1),'fill')
hatch(h(2),'filly')
hatch(h(3),'rare')
hatch(h(4),'Rare')
hatch(h(5),'rarethick')
hatch(h(6),'cross')
hatch(h(7),'plus')
string=char({'fill','filly','rare','Rare','rarethick','cross','plus'});
text(x,y,string);