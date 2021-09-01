function [x,y] = parCircle(theta, c, r)
%PARELLIPSE Returns x and y coordinates of circle parameterized by theta
%   Given the the coordinates of the center c of a circle with radius r,
%   this function computes the x and y coordinates of the points on the
%   circle parameterized by the angles provided in theta, relative to the
%   center of the circle.
%%%
% Compute the coordinates of the circle relative to the center
x = r*cos(theta) + c(1);
y = r*sin(theta) + c(2);
end