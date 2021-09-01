function [x,y] = parEllipse(theta, focus1, focus2, d)
%PARELLIPSE Returns x and y coordinates of ellipse parameterized by theta
%   Given the the coordinates of the foci of an ellipse with major diameter
%   d, this function computes the x and y coordinates of the points on the
%   ellipse parameterized by the angles provided in theta, relative to the
%   center of the ellipse, and aligned with the major axis, so that focus1
%   is at an angle of theta = pi and focus2 is at an angle theta = 0.
%%%
% Compute the center of the circle.
center = 0.5*(focus2 + focus1);
%%%
% Compute the major axis direction and reference angle
diffFoci = focus2 - focus1;
angleMajorDir = atan2(diffFoci(2),diffFoci(1));
%%%
% Get the semi-major and semi-minor lengths
semiMajorLen = d/2;
semiMinorLen = 0.5*sqrt(d^2 - norm(diffFoci)^2);
%%%
% Compute the coordinates of the ellipse relative to the center
xp = semiMajorLen*cos(theta);
yp = semiMinorLen*sin(theta);
%%%
% Transform the ellipse into the proper location and allignment
x = xp*cos(angleMajorDir)-yp*sin(angleMajorDir) + center(1);
y = xp*sin(angleMajorDir)+yp*cos(angleMajorDir) + center(2);
end