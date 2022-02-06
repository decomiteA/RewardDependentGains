function [ boolean_vector ] = BooleanReachedM2Bis( input_structure )
%BooleanReached is a function that determines whether the target has been
%reached or not during a given trial
%Antoine DE COMITE

boolean_vector = zeros(length(input_structure.vector_TP),1);

for jj = 1 : length(boolean_vector)
    if ~(cellfun(@isempty,input_structure.array_timing(jj,5)))
        strint = input_structure.array_timing{jj,5};
        if (strint(1)=='f')
            str = input_structure.array_timing{jj,8};
            boolean_vector(jj) = (str(2)=='a');
        end
    end
end

end

