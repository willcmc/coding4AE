
%%%%%%% 1 and 2
% dx = pi/10;
% x = -2*pi:dx:2*pi;
% f = sin(x);

% df = mynumdiff(f, x, 'backward');
% integral = mynumintegrate(f, x, 'simpson');

%%%%%%% 3
% h = 0.001
% dom = [[-1 1]; [-1,1]];
% Xlist = dom(1,1):h:dom(1,2);
% Ylist = dom(2,1):h:dom(2,2);

% [X,Y] = meshgrid(Xlist, Ylist);

% c = 0.01;
% U = c*X./(X.^2 + Y.^2);
% V = c*Y./(X.^2 + Y.^2);

% knots = [0 0; 0.5, 0.7]

% angles = linspace(0, 2*pi, 100);
% radius = 0.5;
% xCenter = 0;
% yCenter = 0;

% x0 = radius * cos(angles) + xCenter; 
% y0 = radius * sin(angles) + yCenter;

% [~, n] = size(x0)

% for i = 1:n
%     [strX, strY] = mystreamline(U, V, X, Y, x0(i), y0(i));
%     plot(strX, strY, 'b')
%     hold on
% end

% [Xp, Yp] = mypathline(U, V, X, Y, 4, 0.701255, 0.759492);
% plot(Xp, Yp, '*r')


%%%%%%% 4
a = 19;
b = 1+9+1+3+7;
c = 37;
d = -1;

h = 0.05;

dom = [[-1, 1]; [-1,1]; [-1, 1]];

Xlist = dom(1,1):h:dom(1,2);
Ylist = dom(2,1):h:dom(2,2);
Zlist = dom(3,1):h:dom(3,2);

[X, Y, Z] = meshgrid(Xlist, Ylist, Zlist);

for t = 1:100
    u = a*(Y.^2 + Z.^2) + (b + c*t).*X.^3.*Y.*Z + c*exp(d.*X.*t);
    v = a*(Z.^2 + X.^2) + (b + c*t).*X.*Y.^3.*Z + c*exp(d.*Y.*t);
    w = a*(X.^2 + Y.^2) + (b + c*t).*X.*Y.*Z.^3 + c*exp(d.*Z.*t);


    Div = VectCalc.divergence(u, v, w, X, Y, Z);

    z1 = 5;
    contour(X(:,:,z1), Y(:,:,z1), Div(:,:,z1)); hold on

    divMatlab = divergence(X, Y, Z, u, v, w);
    contour(X(:,:,z1), Y(:,:,z1), divMatlab(:,:,z1));

    % title(sprintf('t = %.1f', t));
    % quiver3(X, Y, Z, u, v, w)
    % hold on
    % pause(0.1);
    break;
end