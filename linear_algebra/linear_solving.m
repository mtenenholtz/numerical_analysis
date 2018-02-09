% For more, check out https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution

% Implementation of the forward substitution algorithm for solving a linear system, Ax = B, for x, where A is a lower
% triangular matrix.
function x = forwardsub(A, b)
    n = size(A,1);
    x = zeros(n,1);

    sum = 0;
    for i = 1:n
        sum = (A(i,1:i-1)) * x(1:i-1);
        x(i) = (b(i) - sum) / A(i,i);
    end
end

% Implementation of the back substitution algorithm for solving a linear system, Ax = B, for x, where A is an upper
% triangular matrix.
function x = backsub(A, b)
    n = size(A,1);
    x = zeros(n,1);

    sum = 0;
    for i = n:-1:1
        sum = (A(i,i+1:n)) * x(i+1:n);
        x(i) = (b(i) - sum) / A(i,i);
    end
end
