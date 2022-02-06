function [output_vector] = mappingMod(input_vector)
% antoine de Comite - 20 March 2020
output_vector = zeros(length(input_vector),1);
for ii = 1 : length(output_vector)
   switch input_vector(ii)
       case 4
           output_vector(ii) = 1;
       case 5
           output_vector(ii) = 2;
       case 6
           output_vector(ii) = 3;
       case 7
           output_vector(ii) = 1;
       case 8
           output_vector(ii) = 2;
       case 9
           output_vector(ii) = 3;
   end
end


end