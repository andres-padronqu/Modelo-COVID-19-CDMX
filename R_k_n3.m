% Definir las variables iniciales fuera del alcance de F
I0 = 1/1000000;
S0 = 1 - I0;
E0 = 0;
L0 = 0;
G0 = 0;
H0 = 0;
ICU0 = 0;
REC0 = 0;

% Llamar a la función R0 con las variables iniciales como argumento
Y0 = [S0, E0, I0, L0, G0, H0, ICU0, REC0];

% Llamar a la función principal para obtener los resultados
[t, YY] = R_K_n3(0, 100, Y0);

% Graficar la proporción de contagios activos I(t)
figure;
plot(t, YY(:, 3), 'LineWidth', 2);
xlabel('Tiempo (días)');
ylabel('Proporción de contagios activos');
title('Dinámica de contagios activos de COVID-19 en CDMX');
grid on;

% Encontrar el punto máximo de la serie de contagios activos
[max_I, idx_max] = max(YY(:, 3));

% Encontrar la fecha correspondiente al valor máximo
fecha_pico = datetime('28-Feb-2020') + days(t(idx_max));

disp(['La fecha predicha para el pico de la pandemia es: ', datestr(fecha_pico)]);


% Función principal para simular la dinámica de contagios
function [t,YY] = R_K_n3(t0, T, Y0)
    m = 10000;  % Número de pasos de Runge-Kutta
    h = (T-t0)/m;
    t = linspace(t0, T, m+1);
    n = length(Y0);
    YY = zeros(m+1, n);
    YY(1, :) = Y0;

    for i = 1:m
       k1 = F(t(i), YY(i, :));
       k2 = F(t(i) + h/2, YY(i, :) + h*k1/2);
       k3 = F(t(i) + h/2, YY(i, :) + h*k2/2);
       k4 = F(t(i) + h, YY(i, :) + h*k3);
       YY(i+1, :) = YY(i, :) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end

% Función que define el modelo COVID-19 CDMX (F en el método de Runge-Kutta)
function W = F(t, Y)
    % Definir las constantes del modelo
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

    % Definir las ecuaciones del modelo
    W(1) = -(R0/Dinf)*Y(3)*Y(1);
    W(2) = (R0/Dinf)*Y(3)*Y(1) - Y(2)/Dinc;
    W(3) = Y(2)/Dinc - Y(3)/Dinf;
    W(4) = (1 - Pgrave)/Dinf*Y(3) - Y(4)/Dr;
    W(5) = (Pgrave)/Dinf*Y(3) - Y(5)/Dhosp;
    W(6) = Y(5)/Dhosp - (1 - Picu)/Drh*Y(6) - Picu/Dicu*Y(6);
    W(7) = Picu/Dicu*Y(6) - (1 - Pm)/Dricu*Y(7) - Pm/Dm*Y(7);
    W(8) = Y(4)/Dr + (1 - Picu)/Drh*Y(6) + (1 - Pm)/Dricu*Y(7);
end