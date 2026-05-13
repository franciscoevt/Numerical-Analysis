function [u, c, p, m] = adj_fourier(x, y, p, K)
    x = x(:); y = y(:);
    n_points = length(x);
    m = 2*K + 1;
    omega = pi / p;
    A = zeros(n_points, m);

    for i = 1:n_points
        A(i, 1) = 1.0;
        idx = 2;
        for k = 1:K
            arg = k * omega * x(i);
            A(i, idx) = cos(arg);
            A(i, idx+1) = sin(arg);
            idx = idx + 2;
        end
    end

    % Resolution using internal QR Factorization (minimum norm)
    c = A \ y;
    u = @(t) eval_adj_fourier(t, c, omega, K);
end

function v = eval_adj_fourier(t, c, omega, K)
    t = t(:);
    v = c(1) * ones(size(t));
    idx = 2;
    for k = 1:K
        arg = k * omega * t;
        v = v + c(idx) * cos(arg) + c(idx+1) * sin(arg);
        idx = idx + 2;
    end
end
