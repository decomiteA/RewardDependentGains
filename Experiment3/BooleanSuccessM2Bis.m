function [ boolean_vector,idx_success] = BooleanSuccessM2Bis( input_structure )
%BooleanSuccess is a function that determines whether the trial has been
%successed or not
%
% Antoine De Comite - v1.0
% 2nd of January 2019
boolean_vector = zeros(length(input_structure.vector_TP),1);
length_non_empty = zeros(length(input_structure.vector_TP),1);
for jj = 1 : length(boolean_vector)
    if isempty(find(cellfun(@isempty,input_structure.array_timing(jj,:)),1))
        length_non_empty(jj) = length(input_structure.matrix_timing(jj,:));
    else
    length_non_empty(jj) = min(find(cellfun(@isempty,input_structure.array_timing(jj,:)),1),length(input_structure.matrix_timing(jj,:)));
    strsuccess = input_structure.array_timing{jj,length_non_empty(jj)-2};
    if strsuccess(1) == 's'
        boolean_vector(jj) = 1;
        
    end
end
idx_success = length_non_empty-2;

end

