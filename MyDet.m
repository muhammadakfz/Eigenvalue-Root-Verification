function d = MyDet(A)

    [~, n] = size(A);

    [~, U, ~, tukar] = MyLU(A);

    det_U = 1;
    for i = 1:n
        det_U = det_U * U(i, i);
    end

    d = tukar * det_U;
end
