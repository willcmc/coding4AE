function integral = mynumintegrate(f, x, method)
%mynumintegrate - Perform numerical integration of f over the grid x
%
% Syntax: F = mynumintegrate(f, x, method)
%   options for method: 'trapezoidal' (default), 'simpson'
%
% Performs numerical integration of f over the grid interval x. Either Simpson's rule or the trapezoidal rule may be invoked.
    arguments
        f
        x
        method = 'trapezoidal'
    end

    options = {'trapezoidal', 'simpson'};

    if ~ismember(method, options)
        error("Specify method as either 'trapezoidal', or 'simpson'.", 3, "Unsupported method invoked.");
    end

    n = size(f,2);

    if size(x,2) ~= n
        error("Provide equal-sized function values and input grid.", 2, "Incompatible sizes of f and x");
    end

    if n <= 2 && strcmp(method, 'simpson')
        error("Try again with a finer grid.", 2, "Too few intervals for composite Simpson's rule.");
    end

    integral = 0;

    if strcmp(method, 'trapezoidal')
        for i=2:(n-1)
            integral = integral + 0.5*(f(i-1)+f(i))*(x(i)-x(i-1));
        end

    elseif strcmp(method, 'simpson')
        if ~rem(n, 2)
            n = n-1;
        end

        for i=2:(n/2)
            integral = integral + f(2*i-2) + 4*f(2*i-1) + f(2*i);
        end

        integral = integral*(x(2)-x(1))/3;
    end
end