function ep = VBA_ExceedanceProb (mu, Sigma, form, verbose, Nsamp)

% legacy code
s = warning ('on');
warning ('*** The function `VBA_ExceedanceProb` is now deprecated. Please see `VBA_exceedanceProbability` for an alternative.') 
warning (s);

options = struct ();

try
    form;
catch
    form = 'gaussian';
end

if form(end) == '2'
    form = form(1:end-1);
    options.method = 'analytical';
else
    options.method = 'sampling';
end

try
    options.verbose = verbose;
end

try
    options.nSamples = Nsamp;
end

switch form
    case 'gaussian'
        assert(~isempty(Sigma),'Sigma cannot be empty for Gaussian densities');
        ep = VBA_exceedanceProbability ('Gaussian', mu, Sigma, options);
    case 'dirichlet'
        assert(isempty(Sigma),'Sigma must be empty for Dirichlet densities');
        ep = VBA_exceedanceProbability ('Dirichlet', mu, options);
end



