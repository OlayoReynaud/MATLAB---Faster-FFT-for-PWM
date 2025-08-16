clear; clc; close all;
m_f = [25,50,100,400,1000];

mf_counter = 1;
mean_error_pwmfft = zeros(6,length(m_f));
mean_time_pwmfft = zeros(6,length(m_f));
mean_error_FFT = zeros(6,length(m_f));
mean_time_FFT = zeros(6,length(m_f));

num_samples_pwmfft = 2.^[10,11,12,13,14,15];
num_samples_FFT = 2.^[13,14,15,16,17,18];

for m_f = m_f % Different mf values are evaluated

% num_samples_FFT = 2.^(ceil(log2(m_f)) + [4,5,6,7]);
y = linspace(0,2 * pi * (2^(14+ceil(log2(m_f))) - 1)/(2^(14+ceil(log2(m_f)))),2^(14+ceil(log2(m_f))));
num_levels = 2;
v_c_ref = zeros((num_levels - 1),2^(14+ceil(log2(m_f))));

t_FFT = zeros(1,50 * length(num_samples_FFT));
t_pwmfft = zeros(1,50 * length(num_samples_FFT));

error_FFT = zeros(1,50 * length(num_samples_FFT));
error_pwmfft = zeros(1,50 * length(num_samples_FFT));

for i = 1 : length(num_samples_FFT) % Different number of samples (execution times) are evaluated

K = num_samples_FFT(i);
K_pwmfft = num_samples_pwmfft(i);
cont = 1;

for j = 1 : 30 % 30 points are averaged for each of the above evaluations

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
        order(k) = ceil(rand * m_f/4);
        phi(k) = rand * 2 * pi;
        v_mod_ref = v_mod_ref + A(k) * sin(order(k) * y);
        v_mod = v_mod + A(k) * sin(order(k) * x);
        v_mod_pwmfft = v_mod_pwmfft + A(k) * sin(order(k) * x_pwmfft);
    end

    carrier_phase_array = zeros(1,num_levels - 1);
    for k = 1 : num_levels - 1
        carrier_phase_array(k) = 2 * pi * rand;
    end

    % Reference spectra (FFT with very high sampling):
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

    % tic

    FFT_spectra = 1/length(x) * fft(PWM);

    t_FFT(cont) = toc;

    % Proposed method:

    tic
    pwmfft_spectra = pwmfft(v_mod_pwmfft,m_f,10,num_levels,carrier_phase_array);
    t_pwmfft(cont) = toc;

    % Error:

    error_FFT(cont) = max(abs(abs(FFT_ref_spectra(2 : 4 * m_f)) - abs(FFT_spectra(2 : 4 * m_f))));
    error_pwmfft(cont) = max(abs(abs(FFT_ref_spectra(2 : 4 * m_f)) - abs(pwmfft_spectra(2 : 4 * m_f)/2)));

    cont = cont + 1;

end

% s(mf_counter,i) = scatter(mean(error_Olayo),mean(t_Olayo)*1000,'MarkerEdgeColor',colorscale(mf_counter));
mean_error_pwmfft(i,mf_counter) = mean(error_pwmfft);
mean_time_pwmfft(i,mf_counter) = mean(t_pwmfft);
mean_error_FFT(i,mf_counter) = mean(error_FFT);
mean_time_FFT(i,mf_counter) = mean(t_FFT);

end

mf_counter = mf_counter + 1;

end

%%

figure

box on
grid on
hold on

m_f = [25,50,100,400,1000];
colorscale = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628"];

for i = 1 : length(m_f)
    plot(mean_error_FFT(:,i),mean_time_FFT(:,i)*1000,'Color',colorscale(i),'Marker','^','LineWidth',1.1,'MarkerFaceColor',colorscale(i));
end

for i = 1 : length(m_f)
    plot(mean_error_pwmfft(:,i),mean_time_pwmfft(:,i)*1000,'Color',colorscale(i),'Marker','square','LineWidth',1.1,'MarkerFaceColor',colorscale(i),'LineStyle','--');
end

set(gca,'FontName','Times','XScale','log','YScale','log');
xlabel("Error")
ylabel("Execution time [ms]")
xlim([1e-6,2e-2])
ylim([3e-2,3])

% Legend
for i_leg1 = 1:length(m_f) % Dummy plots for color legend
    h_colors(i_leg1) = plot(NaN, NaN, 's', 'Color', ...
        colorscale(i_leg1), 'MarkerFaceColor', colorscale(i_leg1), ...
        'LineWidth', 1.5, 'DisplayName', ...
        sprintf('$m_f = %g$', m_f(i_leg1)));
end
h_blank = plot(NaN, NaN, 'w-','LineWidth', 1.5, 'DisplayName', '');
h_solid = plot(NaN, NaN, 'k-', ...
        'LineWidth', 1.5, 'DisplayName', 'FFT');
h_dashed = plot(NaN, NaN, 'k--', ...
        'LineWidth', 1.5, 'DisplayName', 'New method');
leg = legend([h_colors,h_blank, h_solid, h_dashed], 'Location',...
    'northeast',  'Interpreter', 'latex');
