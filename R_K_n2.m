I0 = 10^-6;
S0 = 1 - I0;
E0 = 0;
L0 = 0;
G0 = 0;
H0 = 0;
ICU0 = 0;
REC0 = 0;
Y0 = [S0, E0, I0, L0, G0, H0, ICU0, REC0];

[t, YY] = R_K_n3(0, 100, Y0);
plot(t, YY(:, 3))

function [t, YY] = R_K_n3(t0, T, Y0)
    m = 10000;  % Número de pasos de Runge-Kutta
    h = (T - t0) / m;
    t = linspace(t0, T, m + 1);
    n = length(Y0);
    YY = zeros(m + 1, n);
    YY(1, :) = Y0;

    for i = 1:m
        k1 = F(t(i), YY(i, :));
        k2 = F(t(i) + h / 2, YY(i, :) + h * k1 / 2);
        k3 = F(t(i) + h / 2, YY(i, :) + h * k2 / 2);
        k4 = F(t(i) + h, YY(i, :) + h * k3);
        YY(i + 1, :) = YY(i, :) + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end

    function W = F(t, Y)
        n = length(Y);
        W = zeros(1, n);

        % Ejemplo Modelo COVID-19 CDMX
        R0 = 2.83;
        Dinf = 2.9;
        Dinc = 5.2;
        Pgrave = 0.138;
        Dr = 14;
        Dhosp = 4;
        Picu = 0.05;
        Drh = 12;
        Pm = 0.03;
        Dicu = 1;
        Dricu = 7;
        Dm = 8;

        W(1) = -(R0 / Dinf) * Y(3) * Y(1);
        W(2) = (R0 / Dinf) * Y(3) * Y(1) - Y(2) / Dinc;
        W(3) = Y(2) / Dinc - Y(3) / Dinf;
        W(4) = (1 - Pgrave) / Dinf * Y(3) - Y(4) / Dr;
        W(5) = (Pgrave / Dinf) * Y(3) - Y(5) / Dhosp;
        W(6) = Y(5) / Dhosp - (1 - Picu) / Drh * Y(6) - (Picu / Dicu) * Y(6);
        W(7) = (Picu / Dicu) * Y(6) - (1 - Pm) / Dricu * Y(7) - (Pm / Dm) * Y(7);
        W(8) = Y(4) / Dr + (1 - Picu) / Drh * Y(6) + (1 - Pm) / Dricu * Y(7);
    end
end

function y = R0(t)
    c1 = 2.8;
    c2 = 0.5;
    c3 = 2.8;
    I1 = [0, 15];
    I2 = [16, 29];

    if t <= I2(1)
        y = c1;
    elseif t >= I2(1) && t <= I2(2)
        y = c2;
    elseif t >= I2(2) + 1
        y = c3;
    elseif t >= I1(2) && t <= I2(1)
        y = c1 + (c2 - c1) * (t - I1(2)) / (I2(1) - I1(2));
    elseif t >= I2(2) && t <= I2(2) + 1
        y = c3 - (c2 - c3) * (t - I2(2) - 1);
    end
end
