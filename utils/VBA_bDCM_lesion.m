function results=VBA_bDCM_lesion(posterior,out)

%% Get normal profile


%% test lesions of each area
 
nx = size(out.options.inF.A,1);
out.options.verbose = 0;
out.options.displayWin = 0;
out.options.multisession.expanded = 1;

out.dim = VBA_check_struct(out.options.dim, ... ,
    'u',size(out.u,1)       ,...
    'n_t',size(out.u,2)      ...
);

[yp,~,~,~,er] = VBA_simulate (...
        out.dim.n_t,...
        out.options.f_fname,...
        out.options.g_fname,...
        posterior.muTheta,...
        posterior.muPhi,...
        out.u,...
        Inf, ...
        inf(size(posterior.a_sigma)), ...
        out.options,...
        posterior.muX0);
   
        rnormal.y=yp-er;
        
        
parfor iRegion = 1:nx
    effect = ones(1,nx);
    effect(iRegion) = 0;
    posterior_l = posterior;
    posterior_l.muTheta = lesion(posterior,out,effect);
    
      [yp,~,~,~,er] = VBA_simulate (...
        out.dim.n_t,...
        out.options.f_fname,...
        out.options.g_fname,...
        posterior_l.muTheta,...
        posterior.muPhi,...
        out.u,...
        Inf, ...
        inf(size(posterior.a_sigma)), ...
        out.options,...
        posterior.muX0);
   
        rlesion(iRegion).y=yp-er;
end
results.lesion = rlesion;
results.normal = rnormal;
results.out=out;
results.posterior=posterior;


end

function muTheta = lesion(posterior,out,effect)

try 
    inF = out.options.inF{1};
catch
    inF = out.options.inF;
end

muTheta = posterior.muTheta;
nu = out.dim.u;
nx = size(inF.A,1);

A = inF.A;
A(A==1) = inF.indA;
A(1:nx+1:nx^2) = 0; %avoid diagonal elements
for i=1:nx
     [~,~,i_out] = find(A(i,:));
     muTheta(i_out) = muTheta(i_out)*effect(i);
end

for j=1:numel(inF.B)
    B = inF.B{j};
    B(B==1) = inF.indB{j};
    for i=1:nx
        [~,~,i_out] = find(B(i,:));
        muTheta(i_out) = muTheta(i_out)*effect(i);
    end
end

C = inF.C;
C(C==1) = inF.indC;
for i=1:nx
     [~,~,i_out] = find(C(i,:));
     muTheta(i_out) = muTheta(i_out)*effect(i);
end

for j=1:nx
    D = inF.D{j};
    D(D==1) = inF.indD{j};
    for i=1:nx
        [~,~,i_out] = find(D(i,:));
        muTheta(i_out) = muTheta(i_out)*effect(i);
    end
end

end


