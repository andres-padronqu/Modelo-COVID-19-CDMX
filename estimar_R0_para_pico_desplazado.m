function estimar_R0_para_pico_desplazado()

    % Constantes originales
    R0_orig = 2.83;
    Dinf = 2.9; Dinc = 5.2;
    
    % Otras constantes del modelo
    Pgrave = 0.138; Dr = 14; Dhosp = 4;
    Picu = 0.05; Drh = 12; Pm = 0.03; Dicu = 1; Dricu = 7; Dm = 8;
    
    % Definimos el tiempo inicial y el vector de condiciones iniciales
    t0 = 0;
    Y0 = [1-1/1000000, 0, 1/1000000, 0, 0, 0, 0, 0];
    
    % Definimos la función para encontrar el pico de la pandemia
    function pico = encontrar_pico(R0)
        % Definimos la función de derivadas F
        function W = F(t,Y)
            n = length(Y);
            W = zeros(1, n);
            % Modelo COVID-19 CDMX. Definimos las constantes:
            W(1) = -(R0 / Dinf) * Y(3) * Y(1);
            W(2) = (R0 / Dinf) * Y(3) * Y(1) - Y(2) / Dinc;
            W(3) = Y(2) / Dinc - Y(3) / Dinf;
            W(4) = (1 - Pgrave) / Dinf * Y(3) - Y(4) / Dr;
            W(5) = (Pgrave) / Dinf * Y(3) - Y(5) / Dhosp;
            W(6) = Y(5) / Dhosp - (1 - Picu) / Drh * Y(6) - Picu / Dicu * Y(6);
            W(7) = Picu / Dicu * Y(6) - (1 - Pm) / Dricu * Y(7) - Pm / Dm * Y(7);
            W(8) = Y(4) / Dr + (1 - Picu) / Drh * Y(6) + (1 - Pm) / Dricu * Y(7);
        end
        
        % Simulamos la dinámica de contagios
        m = 10000;  % Número de pasos de Runge-Kutta
        T = 200;
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

        % Encontramos el día del pico de la pandemia
        [~, idx_pico] = max(YY(:, 3));
        pico = t(idx_pico);
    end

    % Definimos el día deseado para el pico de la pandemia
    dia_deseado = 50;

    % Definimos un rango de valores para R0
    R0_rango = 2.5:0.01:3.5;

    % Iteramos sobre los valores de R0 para encontrar el valor que produce el pico en el día deseado
    for R0 = R0_rango
        pico = encontrar_pico(R0);
        if pico >= dia_deseado
            fprintf('El valor estimado de R0 para que el pico ocurra %d días después es: %.2f\n', dia_deseado, R0);
            break;
        end
    end

end
