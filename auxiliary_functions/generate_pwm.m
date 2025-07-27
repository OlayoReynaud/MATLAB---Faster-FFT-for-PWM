function pwm = generate_pwm(mod ,car, V_dc_up, V_dc_low)
    % V_dc_up and V_dc_low default to 1 and -1 if not passed

    if nargin < 4
        V_dc_low = -1; % Default value for V_dc_low
    end
    if nargin < 3
        V_dc_up = 1; % Default value for V_dc_up
    end
    vect_size = length(mod);
    pwm = zeros(1, vect_size);
    for i = 1:vect_size
        if (mod(i)>=car(i))
            pwm(i) = V_dc_up;
        else
            pwm(i) = V_dc_low;
        end
    end
end
