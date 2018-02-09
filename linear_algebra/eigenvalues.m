% Requires lu_factorization.m.
% These functions utilize the Rayleigh Quotient to convert eigenvectors to their corresponding eigenvalue.
% https://en.wikipedia.org/wiki/Rayleigh_quotient

% Implementation of the inverse power method algorithm for approximating either the desired eigenvalue or eigenvector.
% Takes in a matrix A to be solved for eigenvalues, and an initial guess lambda0 for the desired eigenvalue. Read more
% here: https://en.wikipedia.org/wiki/Inverse_iteration
function [lambda] = inv_pow_method(A, lambda0)
    x = [1;1;1;1];
    lambda = lambda0;

    it = 1;
    [L,U] = LUFact(A-lambda0*eye(4));
    for i = 1:10
        xt = linsolve(L*U,x);
        x = xt/sqrt(xt'*xt);
        lambda = x'*(A*x);
    end
end


% Implementation of Rayleigh Quotient iteration. Takes in a matrix A and a close-enough initial guess lambda0 to the
% desired eigenvalue, and returns the associated eigenvalue lambda and eigenvector v
function [lambda, v] = rayquot(A, lambda0)
    N = size(A,1);
    v = ones(N,1);

    lambda = l_0;

    for t = 1:10
        v = (A - lambda*eye(N))\v;
        v = v / norm(v,2)
        lambda = v'*(A*v)
    end
end


% Requires gram_schmidt.m
% Implementation of the simultaneous power iteration algorithm, similar to the built-in eig() command in MATLAB. Takes
% in a matrix A and returns all of its associated eigenvectors evects and the eigenvalues evals associated with them.
function [evects, evals] = simultaneous_power_iteration(A)
    N = size(A,2);
    x = eye(N);

    [Q,R] = mgs(A*x);
    for k = 0:10
        [Q,R] = mgs(A*Q);
    end

    evals = R;
    evects = Q;
end


% The next two functions are part of a process for solving for the smallest (i.e. largest negative) eigenvalue, and the
% associated eigenvector, of an input matrix A.
function [lambda, x] = pow_iteration(A)
    x = [1;1;1;1];
    lambda = 0;

    for i = 1:10
        xt = A*x;
        x = xt/norm(xt,1);
        lambda = (x'*(A*x))/(x'*x);
        disp(sprintf("%10d %15.10f",it,lambda));
        it = it + 1;
    end
    if lambda > 0
        [lambda, ~] = alt_inv_pow_method(A,-lambda,actual_smallest);
    end
end

function [lambda, x] = alt_inv_pow_method(A, lambda0, actual)
    x = [1;1;1;1];
    lambda = lambda0;

    [L,U] = LUFact(A-lambda0*eye(4));
    for i = 1:10
        xt = linsolve(L*U,x);
        x = xt/sqrt(xt'*xt);
        lambda = x'*(A*x);
    end
end

