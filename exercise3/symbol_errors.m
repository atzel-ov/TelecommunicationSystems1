function num_of_symbol_errors = symbol_errors(est_X, X)

    num_of_symbol_errors = 0;
    L = length(X);

    for i = 1:L
        if(est_X(i,1) ~= X(i,1) || est_X(i,2) ~= X(i,2))
            num_of_symbol_errors = num_of_symbol_errors + 1;
        end
        
    end

end