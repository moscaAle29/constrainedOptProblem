function xhat = box_projection(x, feasible_set)
% 
% function [xhat] = box_projection(x, mins, maxs)
%
% Function that performs the projection of a vector on the boundaries of an
% n-dimensional box, if the vector is not inside it.
%
% INPUTS:
% x = n-dimensional vector;
% mins = n-dimensional vector where the i-th element is the left boundary
% of the i-th interval chararacterizing the box;
% maxs = n-dimensional vector where the i-th element is the right boundary
% of the i-th interval chararacterizing the box;
%
% OUTPUTS:
% xhat = it is x if x is in the box, otherwise it is the projection
% of x on the boundary of the box.
%
mins = feasible_set(:,1);
maxs = feasible_set(:,2);
xhat = max(min(x, maxs), mins);


end

