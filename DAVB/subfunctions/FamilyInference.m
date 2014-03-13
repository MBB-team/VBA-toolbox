function [p1y,p2y] = FamilyInference(F,ac)

try
    ac;
catch
    ac = 'A';
end

if isequal(ac,'A')
    ind1 = [1,3]; % models containing A(2,1)
    ind2 = [2,3]; % models containing A(1,2)
    all = 1:5;
elseif isequal(ac,'C')
    F = F';
    ind1 = [1,3]; % models containing C(1)
    ind2 = [2,3]; % models containing C(2)
    all = 1:4;
end

mF = max(F(:));

% first model
F1 = F(ind1,:)-mF;
Fno1 = F(setdiff(all,ind1),:)-mF;
p1 = 1/(length(ind1)*size(F,2));
pno1 = 1/(length(setdiff(all,ind1))*size(F,2));
p1y = exp(F1)*p1;
pno1y = exp(Fno1)*pno1;
sump = sum(p1y(:)) + sum(pno1y(:));
p1y = sum(p1y(:))/sump;


% first model
F2 = F(ind2,:)-mF;
Fno2 = F(setdiff(all,ind2),:)-mF;
p2 = 1/(length(ind2)*size(F,2));
pno2 = 1/(length(setdiff(all,ind2))*size(F,2));
p2y = exp(F2)*p2;
pno2y = exp(Fno2)*pno2;
sump = sum(p2y(:)) + sum(pno2y(:));
p2y = sum(p2y(:))/sump;