function [L, U, P, tukar] = MyLU(A)

    [n, ~] = size(A);

    U = A;
    L = eye(n);
    P = eye(n);
    tukar = 1;

    for i = 1:n-1

        [~, idx] = max(abs(U(i:n, i)));
        p = idx + i - 1;

        if p ~= i
            U([i, p], :) = U([p, i], :);
            P([i, p], :) = P([p, i], :);
            tukar = tukar * -1;

            if i > 1
                L([i, p], 1:i-1) = L([p, i], 1:i-1);
            end
        end

        for j = i+1:n
            L(j, i) = U(j, i) / U(i, i);
            U(j, :) = U(j, :) - L(j, i) * U(i, :);
        end
    end
end
