function R0_estimado = estimar_R0_para_pico_retrasado()

    % Parámetros del problema
    dias_retraso_pico = 50; % Cuántos días queremos que se dé el primer pico desde el inicio
    
    % Condiciones iniciales
    I0 = 1/1000000;
    S0 = 1 - I0;
    E0 = 0;
    L0 = 0;
    G0 = 0;
    H0 = 0;
    ICU0 = 0;
    REC0 = 0;
    Y0 = [S0, E0, I0, L0, G0, H0, ICU0, REC0];
    
    % Función para estimar el valor de R0
    function error = calcular_error(R0) % Calcular el error entre el tiempo del pico estimado y el tiempo deseado, dado un valor de R0
        % Definir el modelo con el nuevo valor de R0
        function W = F(t, Y, R0)
            % Constantes del modelo
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
            
            % Definir las ecuaciones del modelo con el nuevo valor de R0
            W(1) = -(R0/Dinf)*Y(3)*Y(1);
            W(2) = (R0/Dinf)*Y(3)*Y(1) - Y(2)/Dinc;
            W(3) = Y(2)/Dinc - Y(3)/Dinf;
            W(4) = (1 - Pgrave)/Dinf*Y(3) - Y(4)/Dr;
            W(5) = Pgrave/Dinf*Y(3) - Y(5)/Dhosp;
            W(6) = Y(5)/Dhosp - (1 - Picu)/Drh*Y(6) - Picu/Dicu*Y(6);
            W(7) = Picu/Dicu*Y(6) - (1 - Pm)/Dricu*Y(7) - Pm/Dm*Y(7);
            W(8) = Y(4)/Dr + (1 - Picu)/Drh*Y(6) + (1 - Pm)/Dricu*Y(7);
        end
        
        % Resolver el modelo con el nuevo valor de R0
        [~, YY] = R_K_n(0, dias_retraso_pico, Y0, @(t, Y) F(t, Y, R0));
        
        % Encontrar el valor máximo de la proporción de contagios activos
        [~, idx_max] = max(YY(:, 3));
        
        % Calcular el error como la diferencia entre el tiempo del pico
        % estimado y el tiempo deseado
        error = dias_retraso_pico - idx_max;
    end

    % Estimar el valor de R0 utilizando búsqueda iterativa
    R0_estimado = fzero(@calcular_error, 2.83); % Valor inicial de R0: 2.83, minimiza el error
end

% Función de Runge-Kutta de cuarto orden
function [t, YY] = R_K_n(t0, T, Y0, F)
    m = 10000;  % Número de pasos de Runge-Kutta
    h = (T - t0) / m;
    t = linspace(t0, T, m + 1);
    n = length(Y0);
    YY = zeros(m + 1, n);
    YY(1, :) = Y0;

    for i = 1:m
        k1 = F(t(i), YY(i, :));
        k2 = F(t(i) + h/2, YY(i, :) + h*k1/2);
        k3 = F(t(i) + h/2, YY(i, :) + h*k2/2);
        k4 = F(t(i) + h, YY(i, :) + h*k3);
        YY(i + 1, :) = YY(i, :) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end
