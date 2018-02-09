% Implementation of the composite trapezoidal rule for integral approximation. Implemented on an example function.
% https://en.wikipedia.org/wiki/Trapezoidal_rule
function In = comptrap(a,b,n)
    h = (b-a)/n;
    sum = 0;

    for i = 1:n-1
        x = a + i*h;
        sum = sum + f(x);
    end

    In = h * (f(a) + 2*sum + f(b))/2;

    function y = f(x)
        y = sin(x)/x;
    end
end


% Implementation of the composite midpoint rule for integral approximation. Implemented on an example function.
% https://en.wikipedia.org/wiki/Riemann_sum
function In = compmid(a,b,n)
    h = (b-a)/n;
    sum = 0;

    for i = 0:n-1
        sum = sum + f((a + h/2)+i*h);
    end

    In = h*sum;

    function y = f(x)
        y = sin(x)/x;
    end
end


% Implementation of the composite Simpson's rule for integral approximation. Implemented on an example function.
% https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_3/8_rule
function In = compsimp(a,b,n)
    h = (b-a)/n;
    x0 = f(a) + f(b);
    x_1 = 0;
    x_2 = 0;

    for i = 1:n-1
        x = a + i*h;
        if mod(i,2) == 1
            x_1 = x_1 + f(x);
        else
            x_2 = x_2 + f(x);
        end
    end

    In = h*(x0 + 2*x_2 + 4*x_1) / 3;

    function y = f(x)
        y = sin(x)/x;
    end
end
