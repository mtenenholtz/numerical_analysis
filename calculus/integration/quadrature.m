% Implementation of two-point Gaussian Quadrature for integral approximation. Implemented on an example function.
% https://en.wikipedia.org/wiki/Gaussian_quadrature
function In = gauss2(a,b,n)
    N = n^2;
    x = linspace(a,b,N);
    h = x(2) - x(1);

    c_1 = 1/2 * h * (1 - 1/sqrt(3));
    c_2 = 1/2 * h * (1 + 1/sqrt(3));
    sum = 0;
    for i = 1:N-1
        sum = sum + f(x(i) + c_1) + f(x(i) + c_2);
    end

    In = 1/2 * h * sum;

    function y = f(x)
        y = exp(-x);
    end
end


% Implementation of three-point Gaussian Quadrature for integral approximation. Implemented on an example function.
% https://en.wikipedia.org/wiki/Gaussian_quadrature
function In = gauss3(a,b,n)
    h = (b-a)/n;
    x = a:h:b;
    tot = 0;

    for i = 1:n
        [pt,w] = lgwt(3,x(i),x(i+1)); % Library code for weights/pts
        tot = tot + sum(f(pt).*w);
    end

    In = tot;

    function y = f(x)
        y = exp(-x);
    end
end

