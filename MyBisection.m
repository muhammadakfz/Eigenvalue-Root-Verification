function lam = MyBisection(f, a, b, tol, max_iter)
    fa = f(a);
    c_old = NaN;

    for k = 1:max_iter
        c = (a + b) / 2;
        fc = f(c);

        if k > 1 && abs((c - c_old) / max(abs(c), eps)) < tol
            lam = c;
            return;
        end

        if fa * fc < 0
            b = c;
        else
            a = c;
            fa = fc;
        end

        c_old = c;
    end

    lam = c;
end
