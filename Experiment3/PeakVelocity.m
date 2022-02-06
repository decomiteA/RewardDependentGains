function [ output_matrix ] = PeakVelocity( speedmatrix,time_vector )
%PEAKVELOCITY is a function that finds the maximal speed during a trial
%   
% Inputs : -speedmatrix is a matrix that contains the speed profiles that
% has to be analysed
%          - time_vector is a vector that contains the time span along
%          which the speedmatrix is provided
% Outputs : - ouput_matrix is a nx2 matrix whose first column contains the
% index of the maximal found velocity (wrt to time_vector) and the second
% column contains the value of the maximal found velocity 
%
% @Antoine De Comite - v1.0.
% 04th of March 2019

output_matrix = zeros(size(speedmatrix,1),2);

for ii = 1 : size(output_matrix,1)
   [output_matrix(ii,2), output_matrix(ii,1)] = max(abs(speedmatrix(ii,:)));
end

