m_f = 100;

num_samples = floor(logspace(3,5,50));
num_samples_pwmfft = floor(logspace(1.5,4,50));
y = linspace(0,2 * pi * (1e6 - 1) * 1e-6,1e6);
num_levels = [2,3,5,8,15];

t_FFT = zeros(1,10 * length(num_samples));
t_pwmfft = zeros(1,10 * length(num_samples));

error_FFT = zeros(1,10 * length(num_samples));
error_pwmfft = zeros(1,10 * length(num_samples));

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
K_pwmfft = num_samples_pwmfft(i);

for j = 1 : 10

    x = linspace(0,2 * pi * (K - 1)/K,K);
    x_pwmfft = linspace(0,2 * pi * (K - 1)/K,K_pwmfft);

    num_harmonics = ceil(8 * rand);
    A = zeros(1,num_harmonics);
    order = zeros(1,num_harmonics);
    phi = zeros(1,num_harmonics);

    v_mod = zeros(1,K);
    v_mod_ref = zeros(1,length(y));
    v_mod_pwmfft = zeros(1,K_pwmfft);
    for k = 1 : num_harmonics
        A(k) = rand/num_harmonics;
        order(k) = rand * m_f/5;
        phi(k) = rand * 2 * pi;
        v_mod_ref = v_mod_ref + A(k) * sin(order(k) * y);
        v_mod = v_mod + A(k) * sin(order(k) * x);
        v_mod_pwmfft = v_mod_pwmfft + A(k) * sin(order(k) * x_pwmfft);
    end

    carrier_phase_array = zeros(1,num_levels - 1);
    for k = 1 : num_levels - 1
        carrier_phase_array(k) = 2 * pi * rand;
    end

    % Reference spectra:

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

    % Proposed method:

    tic
    pwmfft_spectra = pwmfft(v_mod_pwmfft,m_f,5,num_levels,carrier_phase_array);
    t_pwmfft(cont) = toc;

    % Error calculation

    error_FFT(cont) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f)) - abs(FFT_spectra(2 : 5 * m_f))));
    error_pwmfft(cont) = max(abs(abs(FFT_ref_spectra(2 : 5 * m_f)) - abs(pwmfft_spectra(2 : 5 * m_f)/2)));

    cont = cont + 1;

end

end

    scatter(error_pwmfft,t_pwmfft*1000,'MarkerEdgeColor',colorscale(num_levels_counter))

    num_levels_counter = num_levels_counter + 1;

end

hold off
set(gca,'FontName','Times','XScale','log','YScale','log');
xlabel("Error")
ylabel("Execution time [ms]")
xlim([1e-4,1e-1])
ylim([1e-1,1e1])
legend("2 level PWM","3 level PWM","5 level PWM","8 level PWM","15 level PWM")
