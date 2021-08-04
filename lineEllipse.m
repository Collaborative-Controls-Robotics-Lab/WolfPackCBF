function point = lineEllipse(p, slope, focus1, focus2, r, pos)
% Find points on ellipse intersecting a line (given agent 1 or N)

points = zeros(2);
[A,B,C,D,E,F] = ellipseData(focus1, focus2, r);
if abs(slope(1)) < 10^-7 % Vertical line x = const
    x = p(1);
    points(1,:) = [x x];
    a = C;
    b = B*x + E;
    c = A*x^2 + D*x + F;
    % Quadratic formula
    points(2,:) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
    if pos == 1
        if points(2,1) <= points(2,2)
            point = points(:,1);
        else
            point = points(:,2);
        end
    else
        if points(2,1) <= points(2,2)
            point = points(:,2);
        else
            point = points(:,1);
        end
    %if abs(points(2,1) - p(2)) <= abs(points(2,2) - p(2))
    %    point = points(:,1);
    %else
    %    point = points(:,2);
    end
    %if points(2,1) > points(2,2)
    %    point1 = points(:,1);
    %    point2 = points(:,2);
    %else
    %    point1 = points(:,2);
    %    point2 = points(:,1);
    %end
else
    m = slope(2)/slope(1);
    a = A + m*B + m^2*C;
    b = -m*B*p(1) + B*p(2) - 2*m^2*C*p(1) + 2*m*C*p(2) + D + m*E;
    c = m^2*C*p(1)^2 - 2*m*C*p(1)*p(2) + C*p(2)^2 - m*E*p(1) + E*p(2) + F;
    % Quadratic formula
    points(1,:) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
    points(2,:) = m*points(1,:) - m*p(1) + p(2);
    %distances = sqrt((points(1,:) - p(1)).^2 + (points(2,:) - p(2)).^2);
    %idx = find(distances == min(min(distances)));
    %point = points(:,find(distances == min(min(distances))));
    %point1 = points(:,idx);
    %if idx == 1
    %    point2 = points(:,2);
    %else
    %    point2 = points(:,1);
    %end
    if pos == 1
        if points(2,1) <= points(2,2)
            point = points(:,1);
        else
            point = points(:,2);
        end
    else
        if points(2,1) <= points(2,2)
            point = points(:,2);
        else
            point = points(:,1);
        end
    end
end
end