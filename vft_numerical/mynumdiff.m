function df = mynumdiff(f, x, method)
%mynumdiff - computes the numerical derivative of the function f over the grid x
%
% Syntax: df = mynumdiff(f, x, method='central')
%   options for method: 'central' (default), 'forward', 'backward'
%
% Computes the numerical derivative of f, assuming boundedness and analyticity over x.
% default method is the central difference derivative. other invokable methods are forward and backward differences.

    arguments
        f
        x
        method = 'central'
    end
    
    options = {'central', 'forward', 'backward'};

    if ~ismember(method, options)
        error("Specify method as either 'central', 'forward', or 'backward'.", 3, "Unsupported method invoked.");
    end

    n = max(size(f));

    if max(size(f)) ~= n
        error("Provide equal-sized function values and input grid.", 2, "Incompatible sizes of f and x");
    end
    
    df = zeros(n);
    df(1) = NaN;
    df(n) = NaN;

    if strcmp(method, 'central')
        for i=2:(n-1)
            df(i) = (f(i+1) - f(i-1))/(x(i+1) - x(i-1));
        end
    elseif strcmp(method, 'forward')
        for i=2:(n-1)
            df(i) = (f(i+1) - f(i))/(x(i+1) - x(i));
        end
    elseif strcmp(method, 'backward')
        for i=2:(n-1)
            df(i) = (f(i) - f(i-1))/(x(i) - x(i-1));
        end
    else
        error("It's a trap!");
    end

    df = df(:,1);
end