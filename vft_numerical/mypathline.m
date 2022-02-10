function [X, Y] = mypathline(u, v, x, y, t, x0, y0, dom)
%mypathline - Computes the miserable pathlines from x0 y0 till time t
%
% Syntax: [X, Y] = mypathline(u, v, x, y, t, x0, y0, dom)
%
% Determines the trajectory history of a particle injected at x0 y0 
% till time t aka the pathline over domain dom in a 2D flow field
   
    arguments
        u
        v
        x
        y
        t
        x0 = 0
        y0 = 0
        dom = [-1 1; -1 1]
    end

    X = [x0];
    Y = [y0];

    for tau = 1:t
        [~, ix] = min(abs(x(1,:)-X(end)));
        [~, iy] = min(abs(y(:,1)-Y(end)));
        X = [X X(end)+u(iy,ix)];
        Y = [Y Y(end)+v(iy,ix)];
    end
end