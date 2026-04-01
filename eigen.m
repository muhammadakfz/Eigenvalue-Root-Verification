clc; clear;

function [V,D] = MyEigen(A, tol, max_iter)

    n = length(A);

    diag_A = diag(A)'; 
    total_absolut_baris = sum(abs(A), 2)';
    
    simpangan_maksimal = total_absolut_baris - abs(diag_A);
    bawah = min(diag_A - simpangan_maksimal) - 1; 
    atas = max(diag_A + simpangan_maksimal) + 1;

    step = (atas - bawah) / 500; 
    a = bawah;
    
    fa = f_eigen(A, a);
    dfa = f_prime(A, a);

    eigs = [];

    while a < atas && length(eigs) < n
        b = a + step;
        fb = f_eigen(A, b);
        dfb = f_prime(A, b);
        
        if fa * fb <= 0
            lam = MyBisection(@(x) f_eigen(A, x), a, b, tol, max_iter);
            is_multiple = abs(f_prime(A, lam)) < 1e-3;
            eigs = catat_akar(eigs, lam, is_multiple, n);
        end
        
        if dfa * dfb <= 0
            lam_ext = MyBisection(@(x) f_prime(A, x), a, b, tol, max_iter);
            if abs(f_eigen(A, lam_ext)) < 1e-4
                eigs = catat_akar(eigs, lam_ext, true, n);
            end
        end
        
        a = b; 
        fa = fb; 
        dfa = dfb;
    end

    eigs = [eigs, nan(1, n - length(eigs))];

    V = zeros(n);
    D = diag(round(eigs * 10000) / 10000); 

    for i=1:n
        if ~isnan(eigs(i))
            V(:,i) = eigenvector_gauss(A, eigs(i));
        else
            V(:,i) = NaN;
        end
    end
end

function eigs = catat_akar(eigs, lam, is_multiple, max_n)
    is_new = isempty(eigs) || min(abs(eigs - lam)) > 1e-4;
    
    if is_new && length(eigs) < max_n
        eigs = [eigs, lam];
    end
    
    if is_multiple && length(eigs) < max_n
        eigs = [eigs, lam];
    end
end

function f = f_eigen(A, lam)
    M = A - lam * eye(length(A));
    f = MyDet(M);
end

function df = f_prime(A, lam)
    h = 1e-5;
    df = (f_eigen(A, lam + h) - f_eigen(A, lam)) / h;
end

function v = eigenvector_gauss(A, lam)
    n = length(A);
    M = A - lam * eye(n);

    Ared = M(1:n-1, 1:n-1);
    bred = -M(1:n-1, n);
    x = GaussPivot(Ared, bred);

    v = [x; 1];

    v = v / norm(v);
    v = round(v * 10000) / 10000; 
end