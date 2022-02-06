function [ boolean_vector ] = BooleanSuccessM2Bis( input_structure )
%BooleanSuccess is a function that determines whether the trial has been
%successed or not
% Antoine DE COMITE

boolean_vector = zeros(length(input_structure.vector_TP),1);

for jj = 1 : length(boolean_vector)
    if ~(cellfun(@isempty,input_structure.array_timing(jj,5)))
        str1 = input_structure.array_timing{jj,5};
        if (str1(1)=='f')
            str2 = input_structure.array_timing{jj,8};
            if (str2(2)=='a')
                str3 = input_structure.array_timing{jj,9};
                boolean_vector(jj) = (str3(1)=='s');
            else
                boolean_vector(jj) = 0;
            end
        end
    end
end
    
end

