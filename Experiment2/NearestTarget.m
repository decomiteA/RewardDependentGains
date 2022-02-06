function [target] = NearestTarget(end_positions ,targets_positions)
% The aim of this function is to find the nearest target at the end of the movement
%
% Inputs : end_position is a nx2 vector containing the end-position of the hand of the subject
%          targets_position is a 3x2 vector containing the positions of the potential targets
%
% Outputs : is an integer value (1,2 or 3) referring to the nearest target
% respectively for left, central and right targets
%
% @Antoine DE COMITE 
%

dist_to_targets = zeros(size(end_positions,1),3);

for ii = 1 : size(dist_to_targets,2)
   dist_to_targets(:,ii) = sqrt((end_positions(:,1)-targets_positions(ii,1)).^2+(end_positions(:,2)-targets_positions(ii,2)).^2); 
end

[~,target] = min(dist_to_targets,[],2);

end