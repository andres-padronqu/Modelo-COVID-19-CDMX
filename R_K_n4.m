function [t,YY] = R_K_n4(t0, T, Y0, R0)
% Simular numéricamente la dinámica de contagios usando el método de
% Runge-Kutta

% Parámetros del método de Runge-Kutta
m = 10000;  % Número de pasos de Runge-Kutta
h = (T - t0) / m;
t = linspace(t0, T, m+1);
n = length(Y0);
YY = zeros(m+1, n);
YY(1, :) = Y0;

for i = 1:m
    k1 = F(t(i), YY(i, :), R0);
    k2 = F(t(i) + h/2, YY(i, :) + h*k1/2, R0);
    k3 = F(t(i) + h/2, YY(i, :) + h*k2/2, R0);
    k4 = F(t(i) + h, YY(i, :) + h*k3, R0);
    YY(i+1, :) = YY(i, :) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function W = F(t, Y, R0)
% Modelo COVID-19 CDMX. Definimos las constantes:
Dinf = 2.9; Dinc = 5.2; Pgrave = 0.138; Dr = 14; Dhosp = 4;
Picu = 0.05; Drh = 12; Pm = 0.03; Dicu = 1; Dricu = 7; Dm = 8;

% Definimos las ecuaciones del modelo
W = zeros(1, length(Y));
W(1) = -(R0/Dinf) * Y(3) * Y(1);
W(2) = (R0/Dinf) * Y(3) * Y(1) - Y(2)/Dinc;
W(3) = Y(2)/Dinc - Y(3)/Dinf;
W(4) = (1 - Pgrave)/Dinf * Y(3) - Y(4)/Dr;
W(5) = Pgrave/Dinf * Y(3) - Y(5)/Dhosp;
W(6) = Y(5)/Dhosp - (1 - Picu)/Drh * Y(6) - Picu/Dicu * Y(6);
W(7) = Picu/Dicu * Y(6) - (1 - Pm)/Dricu * Y(7) - Pm/Dm * Y(7);
W(8) = Y(4)/Dr + (1 - Picu)/Drh * Y(6) + (1 - Pm)/Dricu * Y(7);
end

function R0_for_peak = find_R0_for_peak(t0, T, Y0, days_after_start)
% Función para encontrar el valor de R0 que produce el pico de la pandemia
% days_after_start días después del inicio de los contagios.

% Definir la función de error que queremos minimizar
error_function = @(R0) abs(days_to_peak(R0, t0, T, Y0) - days_after_start);

% Realizar una búsqueda de valor mínimo
R0_for_peak = fminsearch(error_function, 2.83); % Suponemos un valor inicial para R0
end

function days = days_to_peak(R0, t0, T, Y0)
% Función para calcular los días transcurridos hasta el pico de la pandemia
% dado un valor de R0.

% Definir el número de días en que se evaluará el modelo
num_days = 500;

% Ejecutar el modelo con el valor de R0 dado
[~, YY] = R_K_n(t0, T + num_days, Y0, R0);

% Encontrar el día del pico de la pandemia
[~, idx_peak] = max(YY(:, 3));

% Calcular los días transcurridos hasta el pico de la pandemia
days = idx_peak - 1; % Restamos 1 porque idx_peak incluye el primer día
end

% Parámetros iniciales y ejecución del modelo
% Definir los parámetros iniciales
t0 = 0;
T = 100;
I0 = 1/1000000;
S0 = 1 - I0;
E0 = 0;
L0 = 0;
G0 = 0;
H0 = 0;
ICU0 = 0;
REC0 = 0;

Y0 = [S0, E0, I0, L0, G0, H0, ICU0, REC0];

% Encontrar el valor de R0 que produce el pico de la pandemia 50 días después del inicio de los contagios
R0_for_peak = find_R0_for_peak(t0, T, Y0, 50);

% Ejecutar el modelo con el valor de R0 encontrado
[t, YY] = R_K_n(t0, T, Y0, R0_for_peak);

% Graficar la proporción de contagios activos I(t)
figure;
plot(t, YY(:, 3), 'LineWidth', 2);
xlabel('Tiempo (días)');
ylabel('Proporción de contagios activos');
title('Dinámica de contagios activos de COVID-19 en CDMX');
grid on;

% Encontrar la fecha del pico de la pandemia
[~, idx_peak] = max(YY(:, 3));
fecha_pico = datetime('28-Feb-2020') + days(t(idx_peak));
disp(['La fecha predicha para el pico de la pandemia es: ', datestr(fecha_pico)]);
end
