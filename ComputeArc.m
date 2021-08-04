function arc = ComputeArc(pos1, pos2, radius)
% Computes arclength between two points (column vectors) on a circle with
% given radius
    d = sqrt((pos1(1) - pos2(1))^2 + (pos1(2) - pos2(2))^2);
    th = acos(1 - d^2/(2*radius^2));
    arc = radius*th;
end