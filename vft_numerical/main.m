dx = pi/10;
x = -2*pi:dx:2*pi;
f = sin(x);
% df = mynumdiff(f, x, 'backward');
% integral = mynumintegrate(f, x, 'simpson');

h = 0.001
dom = [[-1 1]; [-1,1]];
Xlist = dom(1,1):h:dom(1,2);
Ylist = dom(2,1):h:dom(2,2);

[X,Y] = meshgrid(Xlist, Ylist);

c = 0.01;
U = c*X./(X.^2 + Y.^2);
V = c*Y./(X.^2 + Y.^2);

knots = [0 0; 0.5, 0.7]

angles = linspace(0, 2*pi, 100);
radius = 0.5;
xCenter = 0;
yCenter = 0;

x0 = radius * cos(angles) + xCenter; 
y0 = radius * sin(angles) + yCenter;

[~, n] = size(x0)

for i = 1:n
    [strX, strY] = mystreamline(U, V, X, Y, x0(i), y0(i));
    plot(strX, strY, 'b')
    hold on
end

[Xp, Yp] = mypathline(U, V, X, Y, 4, 0.701255, 0.759492);
plot(Xp, Yp, '*r')
