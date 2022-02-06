function [ boolean_vector, idx_reached ] = BooleanReachedM2Bis( input_structure )
%BooleanReached is a function that determines whether the target has been
%reached or not during a given trial
%
% Antoine De Comite - v1.0
% 2nd of January

boolean_vector = zeros(length(input_structure.vector_TP),1);
length_non_empty = zeros(length(input_structure.vector_TP),1);
for jj = 1 : length(boolean_vector)
    if isempty(find(cellfun(@isempty,input_structure.array_timing(jj,:)),1))
        length_non_empty(jj) = length(input_structure.matrix_timing(jj,:));
    else
        length_non_empty(jj) = min(find(cellfun(@isempty,input_structure.array_timing(jj,:)),1),length(input_structure.matrix_timing(jj,:)));
        strreached = input_structure.array_timing{jj,length_non_empty(jj)-3};
        if strreached(2) == 'a'
            boolean_vector(jj) = 1;
        end
    end
end
idx_reached = length_non_empty-3;
end

