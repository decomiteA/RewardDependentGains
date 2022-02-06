function [angle] = initialAngle(TotKine,TimeVector)
%   Antoine De Comite 
% 26th February 2020
% The aim of this function is to determine the initial angle taken by the
% hand during movement.


angle = zeros(size(TotKine,1),1);
for ii = 1 : length(angle)
   angle(ii) = atan((TotKine(ii,TimeVector==0,1)-TotKine(ii,TimeVector==-100,1))/(TotKine(ii,TimeVector==0,2)-TotKine(ii,TimeVector==-100,2))); 
end

end

