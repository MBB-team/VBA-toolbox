function options = checkSources(options)

    if ~isfield(options.sources)
        options.sources(1).type=options.binomial;
    end