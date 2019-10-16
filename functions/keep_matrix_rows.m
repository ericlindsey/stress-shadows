function A = keep_matrix_rows(A,components)
    % Function keeps only some rows of the input matrix
    % if 'components' is a number k in (1,2,3), the matrix is sliced via
    % A = A(k:3:end,:)
    % if 'components' is a 2-element list, the 3rd missing component is dropped.
    % Eric Lindsey, 2018
    
    if length(setdiff([1,2,3],components)) + length(components) ~= 3
        error('Error: components values not recognized! Must be one or more elements of the list [1,2,3].')
    end

    % keep only some of the GPS components
    if length(components) == 2
        % keep 2, drop 1 component.
        % find missing component from the list, drop this:
        droprow = setdiff([1,2,3],components);
        A(droprow:3:end, :) = [];
        
    elseif length(components) == 1
        % keep 1 component only
        A = A(components:3:end, :);
        
    end % else, keep all 3, so we do nothing
    
end