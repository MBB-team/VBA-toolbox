function A = allcomb(sets)



% check for empty inputs
q = ~cellfun('isempty',sets) ;

%
sets = cellfun(@unique,sets,'UniformOutput',false);

ni = sum(q) ;
ii = ni:-1:1 ;


if ni==0,
    A = [] ;
else
    args = sets(q) ;
    if ni==1,
        A = args{1}(:) ;
    else
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(ni+1,A{:}),[],ni) ;
    end
end

