function results=VBA_lesionAnalysis(posterior,out)

%% Get normal profile

[results.normal.kernel, results.normal.y] = find_kernel(posterior,out) ;

%% test lesions of each area
nx = size(out.options.inF.A,1);
out.options.verbose = 0;
out.options.displayWin = 0;
parfor iRegion = 1:nx
    effect = ones(1,nx);
    effect(iRegion) = 0;
    posterior_l = posterior;
    posterior_l.muTheta = lesion(posterior,out,effect);
    
    [rlesion(iRegion).kernel, rlesion(iRegion).y] = find_kernel(posterior_l,out) ;
end
results.lesion = rlesion;
results.out=out;
results.posterior=posterior;


% VBA_lesionAnalysisDisplay(results);

end

function muTheta = lesion(posterior,out,effect)

muTheta = posterior.muTheta;
nu = out.dim.u;
nx = size(out.options.inF.A,1);

A = out.options.inF.A;
A(A==1) = out.options.inF.indA;
A(1:nx+1:nx^2) = 0; %avoid diagonal elements
for i=1:nx
     [~,~,i_out] = find(A(i,:));
     muTheta(i_out) = muTheta(i_out)*effect(i);
end

for j=1:nu
    B = out.options.inF.B{j};
    B(B==1) = out.options.inF.indB{j};
    for i=1:nx
        [~,~,i_out] = find(B(i,:));
        muTheta(i_out) = muTheta(i_out)*effect(i);
    end
end

C = out.options.inF.C;
C(C==1) = out.options.inF.indC;
for i=1:nx
     [~,~,i_out] = find(C(i,:));
     muTheta(i_out) = muTheta(i_out)*effect(i);
end

for j=1:nx
    D = out.options.inF.D{j};
    D(D==1) = out.options.inF.indD{j};
    for i=1:nx
        [~,~,i_out] = find(D(i,:));
        muTheta(i_out) = muTheta(i_out)*effect(i);
    end
end

end


