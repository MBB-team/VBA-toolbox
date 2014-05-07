function CUMF = fcdf(x,nu1,nu2);
%FCDF CDF of F-distribution.
%

if ~isequal(size(x),size(nu1)) 
	if isequal(size(nu1),[1 1])
		nu1 = repmat(nu1,size(x));	
	else
		error('nu1 and x must be same size.');
	end
	
elseif ~isequal(size(x),size(nu2))
	
	if isequal(size(nu1),[1 1])
		nu2 = repmat(nu1,size(x));	
	else
		error('nu1 and x must be same size.');
	end
	
end


CUMF = zeros(size(x));
k = ~CUMF;

ind = find(x <= 0);
CUMF(ind) = nan;
if ~isempty(ind)
	ind2 = find(x(ind) == 0);
	CUMF(ind(ind2))  = 0;
end

k(ind) = 0;

%CUMF = quad8('FPDF',0,x,[],[],nu1,nu2); 

x(k) = nu2(k)./(nu2(k) + nu1(k).*x(k));
CUMF(k) = 1 - betainc(x(k),nu2(k)/2,nu1(k)/2);