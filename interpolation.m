% ==================================================
% INTERPOLACIÓN TRIGONOMÉTRICA (FOURIER)
% Eliminación Gaussiana con pivoteo parcial
% ==================================================

function [u, c, p, M, es_par] = interpolacion_trigonometrica(x, y, p)
    x = x(:);
    y = y(:);
    m = length(x);
    if length(unique(x)) ~= m
        error('Hay abscisas repetidas.');
    end

    p_minimo = max(abs(x));
    if nargin < 3 || isempty(p)
        p = p_minimo;
    else
        if p < p_minimo
            p = p_minimo;
        end
    end

    if any(x < -p) || any(x > p)
        error('Puntos fuera del intervalo.');
    end

    if mod(m, 2) == 1
        M = (m - 1) / 2;
        es_par = false;
    else
        M = (m - 2) / 2;
        es_par = true;
    end

    omega = pi / p;
    A = zeros(m, m);

    for i = 1:m
        xi = x(i);
        A(i, 1) = 1.0;
        idx = 2;
        for k = 1:M
            arg = k * omega * xi;
            A(i, idx) = cos(arg);
            A(i, idx+1) = sin(arg);
            idx = idx + 2;
        end
        if es_par
            A(i, m) = cos((M+1) * omega * xi);
        end
    end

    % Eliminación Gaussiana
    A_gauss = A;
    y_gauss = y;
    n = m;

    for k = 1:n-1
        fila_max = k;
        for i = k+1:n
            if abs(A_gauss(i, k)) > abs(A_gauss(fila_max, k))
                fila_max = i;
            end
        end
        if fila_max ~= k
            A_gauss([k, fila_max], :) = A_gauss([fila_max, k], :);
            y_gauss([k, fila_max]) = y_gauss([fila_max, k]);
        end
        for i = k+1:n
            factor = A_gauss(i, k) / A_gauss(k, k);
            for j = k:n
                A_gauss(i, j) = A_gauss(i, j) - factor * A_gauss(k, j);
            end
            y_gauss(i) = y_gauss(i) - factor * y_gauss(k);
        end
    end

    c = zeros(n, 1);
    for i = n:-1:1
        suma = 0;
        for j = i+1:n
            suma = suma + A_gauss(i, j) * c(j);
        end
        c(i) = (y_gauss(i) - suma) / A_gauss(i, i);
    end

    u = @(x_eval) evaluar(x_eval, c, omega, M, es_par, m);
end

function resultado = evaluar(x_eval, c, omega, M, es_par, m)
    x_eval = x_eval(:);
    resultado = c(1) * ones(size(x_eval));
    idx = 2;
    for k = 1:M
        arg = k * omega * x_eval;
        resultado = resultado + c(idx) * cos(arg);
        resultado = resultado + c(idx+1) * sin(arg);
        idx = idx + 2;
    end
    if es_par
        resultado = resultado + c(m) * cos((M+1) * omega * x_eval);
    end
end
