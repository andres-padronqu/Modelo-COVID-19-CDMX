function estimacion_R0_para_pico()
    % Parámetros iniciales
    t0 = 0; % Día de inicio de los contagios (28 de febrero de 2020)
    t_pico = 50; % Día deseado para el pico de la pandemia

    % Estimación de R0
    R0 = -1 / (t_pico - t0);

    fprintf('El valor estimado de R0 para que el pico ocurra en el día %d es: %.4f\n', t_pico, R0);
end


