function [X, Y] = mystreamline(u, v, x, y, x0, y0, dom)
%mystreamline - Computes streamlines for a given 2D velocity field over specified grid
%
% Syntax: streamlines = mystreamline(u, v, x, y)
%
% Returns the streamlines for a given 2D u v component flow over an x y grid

    arguments
        u
        v
        x
        y
        x0 = 0
        y0 = 0
        dom = [-1 1; -1 1]
    end

    % dy/dx = v/u

    dydx = v./u;

    h = x(1, 2) - x(1,1);

    Xp = [x0];
    Yp = [y0];

    Xn = [x0];
    Yn = [y0];

    while Xp(end) <= dom(1,2) && Xp(end) >= dom(1,1) && Yp(end) <= dom(2,2) && Yp(end) >= dom(2,1)
        [~, ix] = min(abs(x(1,:)-Xp(end)));
        [~, iy] = min(abs(y(:,1)-Yp(end)));
        Xp = [Xp Xp(end)+h];
        Yp = [Yp Yp(end)+h*dydx(iy, ix)];
    end

    while Xn(1) <= dom(1,2) && Xn(1) >= dom(1,1) && Yn(1) <= dom(2,2) && Yn(1) >= dom(2,1)
        [~, ix] = min(abs(x(1,:)-Xn(1)));
        [~, iy] = min(abs(y(:,1)-Yn(1)));
        Xn = [Xn(1)-h Xn];
        Yn = [Yn(1)-h*dydx(iy, ix) Yn];
    end

    X = [Xn Xp];
    Y = [Yn Yp];
end