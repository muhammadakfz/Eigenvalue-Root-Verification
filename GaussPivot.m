function x = GaussPivot(A, b)

n = length(b);

Ab = [A b];

for k = 1:n-1
    [~, idx] = max(abs(Ab(k:n, k)));
    p = idx + k - 1;

    if p ~= k
        Ab([k, p], :) = Ab([p, k], :);
    end

    for i = k+1:n
        factor = Ab(i, k) / Ab(k, k);
        Ab(i, k:n+1) = Ab(i, k:n+1) - factor * Ab(k, k:n+1);
    end
end

x = zeros(n, 1);
for i = n:-1:1
    if i < n
        s = Ab(i, i+1:n) * x(i+1:n);
    else
        s = 0;
    end
    
    x(i) = (Ab(i, n+1) - s) / Ab(i, i);
end
end
