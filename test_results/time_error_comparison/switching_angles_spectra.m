function c = switching_angles_spectra(v_mod,m_f,max_order)

% v_mod : Señal moduladora

% m_f : Índice de modulación en frecuencia

% max_order: Orden del armónico de mayor frecuencia considerado



j = sqrt(-1);
x = linspace(0,2 * pi,length(v_mod));

% Se obtienen los ángulos de conmutación:

v_c = - sawtooth(m_f * x,0.5);
alpha = zeros(1,2 * m_f);
cont = 0;
flag = 0;

for k = 2 : length(x) - 1
    % Bloque de corrección (para ángulos de conmutación no detectados):
    if (v_c(k + 1) < v_c(k)) && (v_c(k - 1) < v_c(k)) || (k == length(x) - 1) % Fin del intervalo de conmutación
        if flag == 3
            flag = 0;
            continue
        end
        if flag == 2
            alpha(cont + 1) = alpha(cont);
            alpha(cont) = x(k + 1) - 2 * pi/m_f;
            cont = cont + 1;
            flag = 0;
            continue
        end

        if flag == 1
            cont = cont + 1;
            alpha(cont) = x(k);
            flag = 0;
            continue
        end

        if flag == 0
            cont = cont + 1;
            alpha(cont) = x(k) - 2  * pi/m_f;
            cont = cont + 1;
            alpha(cont) = x(k);
            flag = 0;
            continue
        end

        % if flag > 0
        %     error("Error en el bloque de corrección.");
        % end

    end

    % Comparación:

    if ((v_c(k + 1) > v_mod(k + 1)) && ((v_c(k) < v_mod(k)) || (v_c(k) == v_mod(k)))) && (flag ~= 3)
        cont = cont + 1;
        alpha(cont) = x(k);
        flag = flag + 2;
        continue
    end

    if ((v_c(k + 1) < v_mod(k + 1)) && ((v_c(k) > v_mod(k)) || (v_c(k) == v_mod(k)))) && (flag ~= 3)
        cont = cont + 1;
        alpha(cont) = x(k);
        flag = flag + 1;
    end
end

for k = 2 : length(x) - 1

    if ((v_c(k + 1) > v_mod(k + 1)) && ((v_c(k) < v_mod(k)) || (v_c(k) == v_mod(k)))) && (flag ~= 3)
        cont = cont + 1;
        alpha(cont) = x(k);
    end

    if ((v_c(k + 1) < v_mod(k + 1)) && ((v_c(k) > v_mod(k)) || (v_c(k) == v_mod(k)))) && (flag ~= 3)
        cont = cont + 1;
        alpha(cont) = x(k);
    end

end

alpha = alpha(1 : cont); % Se reduce el tamaño del vector alpha en caso de haber omitido conmutaciones



% Se obtiene el espectro de la señal PWM a partir de los ángulos de conmutación

c = zeros(1,max_order);

for n = 1 : max_order
    for k = 1 : cont/2
        c(n) = c(n) + exp(- n * alpha(2 * k) * j) - exp(- n * alpha(2 * k - 1) * j);
    end
    c(n) = j/(4 * pi * n) * c(n);
end

end
