clear; clc;;
parentDir = pwd; % Parent directory
p1 = fileparts(parentDir);  % One level up
p2 = fileparts(p1);  % Two level up
addpath(fullfile(p2, 'algorithm_functions/'));

m_f = 1000;

% num_samples = floor(logspace(3,5,100));
% num_samples_Olayo = floor(logspace(1.5,4,100));

num_samples = 2.^(13:17);
num_samples_Olayo = 2.^(5:13);


y = linspace(0,2 * pi,2^18);
v_c_ref = sawtooth(m_f * y,0.5);

numb_sigs = 1e3; % Total number of signals computed

t_FFT = zeros(1,numb_sigs);
t_Olayo = zeros(1,numb_sigs);
% t_sw_ang = zeros(1,numb_sigs);

error_FFT = zeros(1,numb_sigs);
error_Olayo = zeros(1,numb_sigs);
% error_sw_ang = zeros(1,numb_sigs);

for i = 1 : numb_sigs
    
    if rem(i, 150) == 0
        disp(num2str(100 * i / numb_sigs) + " %")
    end
    K = num_samples(randi(length(num_samples))); % Random elemnt from list
    K_Olayo = num_samples_Olayo(randi(length(num_samples_Olayo)));
    

    x = linspace(0,2 * pi,K);
    x_Olayo = linspace(0,2 * pi,K_Olayo);

    num_harmonics = ceil(8 * rand);
    A = zeros(1,num_harmonics);
    order = zeros(1,num_harmonics);
    phi = zeros(1,num_harmonics);

    % Código de Olayo para crear las señales
    v_mod_ref = random_signal(num_harmonics, length(y));
    v_mod_Olayo = random_signal(num_harmonics, length(x_Olayo));
    v_mod = random_signal(num_harmonics, length(x));

    % FFT de referencia:    
    PWM = zeros(1,length(y));    
    for k = 1 : length(y)
        if v_mod_ref(k) > v_c_ref(k)
            PWM(k) = 1;
        end
    end

    FFT_ref_spectra = 1/length(y) * fft(PWM);    
    % FFT:

    tic

    v_c = sawtooth(m_f*x, 0.5);
    PWM = zeros(1,K);    
    for k = 1 : K
        if v_mod(k) > v_c(k)
            PWM(k) = 1;
        end
    end    
    FFT_spectra = 1/length(x) * fft(PWM);    
    t_FFT(i) = toc;

    % Método de Olayo:    
    tic
    Olayo_spectra = get_spectrum_extended(v_mod_Olayo,1,0,m_f,1,5);
    t_Olayo(i) = toc;

    % Método de los ángulos de conmutación:    
    % tic
    % Switching_angles_spectra = switching_angles_spectra(v_mod,m_f,5 * m_f);
    % t_sw_ang(i) = toc;

    % Error cometido
    max_idx = min([5 * m_f + 1, length(FFT_spectra),  length(Olayo_spectra)]);
    error_FFT(i) = max(abs(abs(FFT_ref_spectra(2 : max_idx)) - abs(FFT_spectra(2 : max_idx))));
    error_Olayo(i) = max(abs(abs(FFT_ref_spectra(2 : max_idx)) - abs(Olayo_spectra(2 : max_idx))/2));
    %error_sw_ang(i) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f + 1)) - abs(Switching_angles_spectra)));


end

%% Resultados:

figure

hold on

scatter(error_FFT,t_FFT * 1000)
scatter(error_Olayo,t_Olayo * 1000)
%scatter(error_sw_ang,t_sw_ang * 1000)

hold off
grid on
box on


set(gca,'FontName','Times');
xlabel("Error")
ylabel("Execution time (ms)")

%title("Results for m_f = " + num2str(m_f))

legend('FFT','Proposed method')%,'Switching angles method')


function sig = random_signal(band_width, vect_size)
    % Generates a random singal with amplitude smaller or equal to 1 (it
    % might have harmonics with bigger amplitude).
    sample_length = 2*pi/vect_size;
    y = 0:sample_length:2*pi - sample_length;
    unscaled_signal = zeros(1, vect_size);
    for vm_hrmnc = 1:band_width        
        % Add the harmonics randomly. Scale at the end to have an 
        % amplitude of 1.
        amp = rand(); % Random amplitude
        phase = rand()*2*pi;% Random phase
        unscaled_signal = unscaled_signal + amp*sin(vm_hrmnc*y + phase);
    end
    max_abs = max(abs(unscaled_signal));
    if max_abs == 0
        sig = unscaled_signal; % Zeros
    else
        sig = unscaled_signal/max_abs; % Scale to have range in [-1, 1]
    end    
end
