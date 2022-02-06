function [ y ] = customMean( x, dim, exclude)
%CUSTOMMEAN computes the mean value with some special features
%   This function is quite similar to the mean function of matlab. It adds
%   a special feature that allows to exclude certain entries of the
%   computation. 
%   This special features is the following: If the parameter exclude is non 
%   empty, all the values contained in the exclude vector are excluded from
%   the computation of the mean. In fact they are replaced by NaN so that,
%   by calling the mean function (the one due to MathWorks), they are not
%   taken into account in the computation of the mean.
%
%
%   INPUTS : - x is the entity whose mean wants to be computed, if it's a
%                  vector, the output is the mean value of the vector. If
%                  it's a matrix, by default, it computes the mean along
%                  the first dimension of this matrix. Therefore the output
%                  is a matrix whose dimension is the same as the initial
%                  one but whose first dimension is equal to zero.
%               - dim is the dimension along which the means has to be
%                  computed in the case of matrix. This is an optionnal parameter.
%               - exclude is either a scalar of a vector that has not to be
%                 taken into account in the computation of the mean value of
%                 the x input. This is also an optionnal parameter.
%
%   OUTPUTS : - y is the computed mean value of the input vector or matrix
%
% @Antoine De Comite 
% 10-01-2018 / version 1.0


switch nargin
    
    case 1
        y = mean(x);
    case 2
        if (length(dim)~=1 || dim <1)
            error('The second parameter, the dimension of the computation of the mean has to be an integer greater or equal than 1');
        else
            y = mean(x,dim);
        end
    case 3
        for i=1:length(exclude)
            x(x==exclude(i))=NaN;
        end
        y = mean(x,dim,'omitnan');
end


end

