function response = simulResponses(lambdas)

alter = size(lambdas,1);
n = size(lambdas,2);

% == normalization 
max_lambda = max(sum(lambdas,2)) ;
if max_lambda >1
    lambdas = lambdas / max_lambda ;
end



response = zeros(alter,n-1) ;

for i=1:n-1
    response(:,i) = getRespAt(lambdas,i);
end

ok=false;
while ~ok 
    rs = find(sum(response)) ;
    if numel(rs)==0
        response = simulResponses(lambdas);
    elseif numel(rs)==1
        ok=true;
    else
        test = rs(randi(numel(rs))) ;
        response(:,test) = 1*(response(:,test) & getRespAt(lambdas,test));
    end
    
end

end
function r=getRespAt(lambdas,t)
    l = .5*(lambdas(:,t) + lambdas(:,t+1))';
    r=mnrnd(1,[l 1-sum(l)]) ;
    r=r(1:end-1)';
end


