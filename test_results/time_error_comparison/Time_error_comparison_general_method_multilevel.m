m_f = 100;

num_samples = floor(logspace(3,5,50));
num_samples_Olayo = floor(logspace(1.5,4,50));
y = linspace(0,2 * pi * (1e6 - 1) * 1e-6,1e6);
num_levels = [2,3,5,8,15];

t_FFT = zeros(1,10 * length(num_samples));
t_Olayo = zeros(1,10 * length(num_samples));

error_FFT = zeros(1,10 * length(num_samples));
error_Olayo = zeros(1,10 * length(num_samples));

num_levels_counter = 1;

colorscale = ["#CC0000","#CCCC00","#00CC00","#00CCCC","#0000CC"];

figure

hold on
box on
grid on

for num_levels = num_levels

cont = 1;
v_c_ref = zeros((num_levels - 1),1e6);

for i = 1 : length(num_samples)

K = num_samples(i);
K_Olayo = num_samples_Olayo(i);

for j = 1 : 10

    x = linspace(0,2 * pi * (K - 1)/K,K);
    x_Olayo = linspace(0,2 * pi * (K - 1)/K,K_Olayo);

    num_harmonics = ceil(8 * rand);
    A = zeros(1,num_harmonics);
    order = zeros(1,num_harmonics);
    phi = zeros(1,num_harmonics);

    % Código de Olayo para crear las señales
    % v_mod_ref = random_signal(rand * m_f/5, length(y));
    % v_mod_Olayo = random_signal(rand * m_f/5, length(x_Olayo));
    % v_mod = random_signal(rand * m_f/5, length(x));
    % Código de Jaime para crear la señales
    v_mod = zeros(1,K);
    v_mod_ref = zeros(1,length(y));
    v_mod_Olayo = zeros(1,K_Olayo);
    for k = 1 : num_harmonics
        A(k) = rand/num_harmonics;
        order(k) = rand * m_f/5;
        phi(k) = rand * 2 * pi;
        v_mod_ref = v_mod_ref + A(k) * sin(order(k) * y);
        v_mod = v_mod + A(k) * sin(order(k) * x);
        v_mod_Olayo = v_mod_Olayo + A(k) * sin(order(k) * x_Olayo);
    end

    carrier_phase_array = zeros(1,num_levels - 1);
    for k = 1 : num_levels - 1
        carrier_phase_array(k) = 2 * pi * rand;
    end

    % FFT de referencia:

    if mod(num_levels,2) == 0    
        for k = 1 : num_levels - 1
            v_c_ref(k,:) = (sawtooth(m_f * (y + carrier_phase_array(k)/m_f),0.5) + 2 * (k - num_levels/2))/(num_levels - 1);
        end
    else
        for k = 1 : num_levels - 1
            v_c_ref(k,:) = (sawtooth(m_f * (y + carrier_phase_array(k)/m_f),0.5) + 1 + 2 * (k - ceil(num_levels/2)))/(num_levels - 1);
        end    
    end

    PWM = zeros(1,length(y));    
    for k1 = 1 : length(y)    
        for k2 = 1 : num_levels - 1       
            if v_mod_ref(k1) > v_c_ref(k2,k1)    
               PWM(k1) = PWM(k1) + 1/(num_levels - 1);     
            end
            if v_mod_ref(k1) < v_c_ref(k2,k1)
                PWM(k1) = PWM(k1) - 1/(num_levels - 1);
            end
        end
    
    end

    FFT_ref_spectra = 1/length(y) * fft(PWM);

    % FFT:

    tic

    v_c = sawtooth(m_f * x,0.5);
    PWM = zeros(1,K);

    for k = 1 : K
        if v_mod(k) > v_c(k)
            PWM(k) = 1;
        end
    end

    v_c = zeros((num_levels - 1),K);
    
    if mod(num_levels,2) == 0
    
        for k = 1 : num_levels - 1
            v_c(k,:) = (sawtooth(m_f * (x + carrier_phase_array(k)/m_f),0.5) + 2 * (k - num_levels/2))/(num_levels - 1);
        end    
    
    else
    
        for k = 1 : num_levels - 1
            v_c(k,:) = (sawtooth(m_f * (x + carrier_phase_array(k)/m_f),0.5) + 1 + 2 * (k - ceil(num_levels/2)))/(num_levels - 1);
        end
    
    end

    PWM = zeros(1,K);
    
    for k1 = 1 : K
    
        for k2 = 1 : num_levels - 1       
            if v_mod(k1) > v_c(k2,k1)
    
               PWM(k1) = PWM(k1) + 1/(num_levels - 1); 
    
            end
    
            if v_mod(k1) < v_c(k2,k1)
    
                PWM(k1) = PWM(k1) - 1/(num_levels - 1);
    
            end
        end
    
    end

    FFT_spectra = 1/length(x) * fft(PWM);

    t_FFT(cont) = toc;

    % Método de Olayo:

    tic
    % Olayo_spectra = fft_for_pwm(v_mod_Olayo,1,0,m_f,1,5);
    Olayo_spectra = pwmfft(v_mod_Olayo,m_f,5,num_levels,carrier_phase_array);
    t_Olayo(cont) = toc;

    % Error cometido

    error_FFT(cont) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f)) - abs(FFT_spectra(2 : 5 * m_f))));
    error_Olayo(cont) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f)) - abs(Olayo_spectra(2 : 5 * m_f)/2)));
    % error_sw_ang(cont) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f)) - abs(Switching_angles_spectra)));

    cont = cont + 1;

end

end

    scatter(error_Olayo,t_Olayo*1000,'MarkerEdgeColor',colorscale(num_levels_counter))

    num_levels_counter = num_levels_counter + 1;

end

hold off
set(gca,'FontName','Times','YScale','log');
xlabel("Error")
ylabel("Execution time [ms]")
xlim([0,0.1])
ylim([0,5])
legend("2 level PWM","3 level PWM","5 level PWM","8 level PWM","15 level PWM")

%% Resultados:
% 
% figure
% 
% hold on
% 
% scatter(error_FFT,t_FFT * 1000)
% scatter(error_Olayo,t_Olayo * 1000)
% % scatter(error_sw_ang,t_sw_ang * 1000)
% 
% hold off
% grid on
% box on
% 
% 
% set(gca,'FontName','Times');
% xlabel("Error")
% ylabel("Execution time [ms]")
% xlim([0,0.01])
% ylim([0,4])
% 
% % title("Resultados para m_f = " + num2str(m_f))
% 
% legend('FFT','Proposed method','Método de los ángulos de conmutación')
% 
% 
% function sig = random_signal(band_width, vect_size)
%     % Generates a random singal with amplitude smaller or equal to 1 (it
%     % might have harmonics with bigger amplitude).
%     sample_length = 2*pi/vect_size;
%     y = 0:sample_length:2*pi - sample_length;
%     unscaled_signal = zeros(1, vect_size);
%     for vm_hrmnc = 1:band_width        
%         % Add the harmonics randomly. Scale at the end to have an 
%         % amplitude of 1.
%         amp = rand(); % Random amplitude
%         phase = rand()*2*pi;% Random phase
%         unscaled_signal = unscaled_signal + amp*sin(vm_hrmnc*y + phase);
%     end
%     max_abs = max(abs(unscaled_signal));
%     if max_abs == 0
%         sig = unscaled_signal; % Zeros
%     else
%         sig = unscaled_signal/max_abs; % Scale to have range in [-1, 1]
%     end    
% end
