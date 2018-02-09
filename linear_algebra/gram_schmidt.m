% Note: gram_schmidt() is only for demonstration. Always use modified_gram_schmidt() in practice.
% The Gram-Schmidt process for performing QR factorization. Takes in a matrix A, and returns an orthogonal matrix Q and
% an upper triangular matrix R. Check more out here: https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
function [Q, R] = gram_schmidt(A)
    n = size(A, 2);

    for i = 1:n
        sum = 0;
        for j = 1:i-1
            prod = A(:,i)' * Q(:,j);
            sum = sum + prod*Q(:,j);
            R(j, i) = prod;
        end
        temp = A(:,i) - sum;
        normal = norm(temp, 2);
        Q(:,i) = temp/normal;
        R(i, i) = normal;
    end
end

% The Gram-Schmidt process, modified to make it more numerically stable. You should always use this in practice over
% the original algorithm, as this fixes the floating point errors that end up creating runaway error.
function [Q, R] = modified_gram_schmidt(A)
    n = size(A, 2);

    for i = 1:n
        temp = A(:,i);
        for j = 1:i-1
            prod = temp' * Q(:,j);
            temp = temp - prod * Q(:, j);
            R(j,i) = prod;
        end
        Q(:,i) = temp/norm(temp, 2);
    end
end
