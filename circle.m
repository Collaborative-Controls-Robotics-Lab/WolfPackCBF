function fun = circle(x0,y0,r)
%hold on
%th = 0:pi/50:2*pi;
%xs = r * cos(th) + x0;
%ys = r * sin(th) + y0;
fun = @(x,y) (x - x0).^2 + (y - y0).^2 - r^2;
%hold off
end