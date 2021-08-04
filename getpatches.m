function [xs,ys] = getpatches(x,length)
% Get patches for pretty plotting
% Returns two matrices of x-coordinates and y-coordinates responsively. 
% Each column of the matrix is the coordinates
xs = zeros(size(x,2), 3);
ys = zeros(size(x,2), 3);
for ii = 1:size(x,2)
    u = [cos(x(3,ii)) sin(x(3,ii))];
    uperp = [-sin(x(3,ii)), cos(x(3,ii))];
    xs(ii,:) = [x(1,ii) + length/3*uperp(1) x(1,ii) - length/3*uperp(1) ...
        x(1,ii) + length*u(1)];
    ys(ii,:) = [x(2,ii) + length/3*uperp(2) x(2,ii) - length/3*uperp(2) ...
        x(2,ii) + length*u(2)];
end
xs = xs';
ys = ys';
end