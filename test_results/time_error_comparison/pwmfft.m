function PWM_spectra = pwmfft(v_mod,m_f,num_lobes,num_levels,carrier_phase)

    N = length(v_mod);
    PWM_spectra = zeros(1,num_lobes * m_f + ceil(N/2));
    f = zeros(1,N);

    for m = 0 : num_lobes
        % f calculation:
        if num_levels == 2
            f = exp(- m * carrier_phase * 1i) * 4/(m*pi) * sin(pi*m/2 * (1 + v_mod));
        else        
            for k = 1 : N    
                kappa = 1 + floor(0.9999 * (v_mod(k) + 1) * (num_levels - 1)/2);
                v_mod_kappa = v_mod(k) * (num_levels - 1) + 2 * (ceil((num_levels - 1)/2) - kappa) + mod(num_levels,2);
                if m == 0
                    f(k) = (4 * (kappa - 1) + 2 * (1 + v_mod_kappa))/(num_levels - 1);                
                else    
                    f(k) = exp(- m * carrier_phase(kappa) * 1i) * 4/(m*pi) * sin(pi*m/2 * (1 + v_mod_kappa))/(num_levels - 1);
                end
            end
        end

        F = ifft(f);
        % Harmonics placement:
        if m == 0
            PWM_spectra(1 : floor(N/2)) = F(1 : floor(N/2));
        else
            PWM_spectra(m*m_f + 1) = PWM_spectra(m*m_f + 1) + F(1);            
            if floor(N/2) >= length(PWM_spectra) - m*m_f - 1
                PWM_spectra(m*m_f + 2 : length(PWM_spectra)) = F(2 : length(PWM_spectra) - m*m_f);
            else
                PWM_spectra(m*m_f + 2 : m*m_f + floor(N/2)) = F(2 : floor(N/2));
            end
            if floor(N/2) >= m*m_f
                PWM_spectra(1 : m*m_f) = PWM_spectra(1 : m*m_f) + F(length(F) - m*m_f + 1 : end);
            else
                PWM_spectra(m*m_f - floor(N/2) : m*m_f) = PWM_spectra(m*m_f - floor(N/2) : m*m_f) + F(length(F) - floor(N/2) : end);
            end
        end
    end
    PWM_spectra(1) = PWM_spectra(1) - 2;
end

%% Input parameters:

% v_mod : Modulating signal
% m_f = Frequency index
% num_lobes : Number of harmonic lobes to be calculated
% num_levels : Number of levels of the PWM signal
% carrier_phase : Array of size (num_levels - 1) containing the phase of each carrier signal
