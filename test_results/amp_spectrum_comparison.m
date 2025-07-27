% Used to test-out the results of the new algorithm. The amplidue spectrums
% are compared with respect to the one obtained by FFT.
clear; close all; clc;

parentDir = pwd; % Parent directory
p1 = fileparts(parentDir);  % One level up
addpath(fullfile(p1, 'auxiliary_functions/'));
addpath(fullfile(p1, 'algorithm_functions/'));

% carrier_type can be symetric, leading_edge or trailing_edge (1, 2 or 3)
carrier_type = 1; 

show_time_domain_plots = true;
numb_reg = 2; % Number of regions computed

Vp = 1; % Amplitude of carrier signal
fm = 50; % Frequency of modulating signal

V_dc_up = 1; % Upper value of output square wave
V_dc_low = -1; % Lower value of output square wave
V_dc = V_dc_up - V_dc_low;
mf = 30;

fc = mf*fm; % Frequency of carrier signal 


% Adjust desired sample length
Tm = 1/fm;
vect_size_raw_fft = 2097152;
vect_size = 2048; %1024; %2048; %524288; %8192;

sample_length_raw_fft = Tm/vect_size_raw_fft;
sample_length = Tm/vect_size;
t_raw_fft = 0:sample_length_raw_fft:(Tm - sample_length_raw_fft);
t = 0:sample_length:(Tm - sample_length);

% Create modulating
y = 2*pi*fm*t;
y_raw_fft = 2*pi*fm*t_raw_fft;
vm = cos(y);
vm_raw_fft = cos(y_raw_fft);

% vm = 0.2 + 0.65*cos(y) - 0.35.*cos(2*y) - 0.3.*cos(3*y);

% vm = sawtooth(y, 1/2);
% vm_raw_fft = sawtooth(y_raw_fft, 1/2);

% vm = cos(y) - 0.5.*cos(3*y);
% vm_raw_fft = cos(y_raw_fft) - 0.16.*cos(3*y_raw_fft);

% vm = 0.7*exp(cos(y)) - 1; % Random 2-pi periodic function in y with a range 
% %                           % between -1 and 1
% vm_raw_fft = 0.7*exp(cos(y_raw_fft)) - 1; 

tic
% Create carrier
if carrier_type == 1 % "symetric"
    vc = Vp.*sawtooth(2*pi*fc*t_raw_fft, 1/2);
elseif carrier_type == 2 % "leading_edge"
    vc = Vp.*sawtooth(2*pi*fc*t_raw_fft - pi, 0);
elseif carrier_type ==  3 % "trailing_edge"
    vc = Vp.*sawtooth(2*pi*fc*t_raw_fft - pi, 1);
end

% Create pwm
pwm = generate_pwm(vm_raw_fft, vc, V_dc_up, V_dc_low);
time_to_simulate_pwm = toc

tic
P1 = amp_spectrum_fft(pwm);
total_time_for_fft = toc + time_to_simulate_pwm
Fs_fft = 1/vect_size_raw_fft; % Samplig frequency FFT
freqs_ffts = Fs_fft*(0:(vect_size_raw_fft/2))/vect_size_raw_fft;


tic
% The number of regions to be computed can be added at the end of the
% following function (as a possitive integer).
pwm_hmncs = get_spectrum_extended(vm, V_dc_up, V_dc_low, mf, carrier_type, numb_reg);
time_for_new_algorithm = toc
output_numb_samples = length(pwm_hmncs); % Should also be vect_size/2

pwm_hmncs(1) = pwm_hmncs(1)/2;
amp_err = P1(1:output_numb_samples) - abs(pwm_hmncs);

% Plots
if show_time_domain_plots
    subplot(2,2,1)
    plot_mod_and_carr(vm_raw_fft, vc, t_raw_fft)
    subplot(2,2,3);
    plot_pwm(t_raw_fft, pwm)
end

if show_time_domain_plots
    subplot(2,2,2)
else
    subplot(3,1,1)
end

plot_amp_spectrum(P1(1:output_numb_samples), 1e-3)
title("Espectro mediante FFT")
xlim([1, output_numb_samples]);

if show_time_domain_plots
    subplot(2,2,4)
else
    subplot(3,1,2)
end

plot_amp_spectrum(pwm_hmncs, 1e-3)
title('Espectro mediante FFT simplificada');
xlim([1, output_numb_samples]);

if not(show_time_domain_plots)
    subplot(3,1,3)
    plot_amp_spectrum(amp_err, 1e-5)
    title("Diferencia entre ambos");
    xlim([1, output_numb_samples]);
end

function plot_mod_and_carr(vm, vc, t)
    Vp = max(vc);
    yyaxis right
    plot(t, vc, 'linewidth', 0.75);
    ylim([-Vp - 0.5*Vp, Vp + 0.5*Vp])
    ylabel("Portadora")

    % Horizontal lines for amplitude of carrier
    hold on
    plot([0, t(end)], [Vp, Vp], '--', 'linewidth', 0.5);
    plot([0, t(end)], [-Vp, -Vp], '--', 'linewidth', 0.5);
%     set(gca, 'ytick', [-Vp, Vp],'yticklabel', ["-Vp"; "Vp"])

    yyaxis left
    plot(t, vm, 'linewidth', 1.5);
    ylabel("Moduladora")

    xlabel('Tiempo');
    ylim([-Vp - 0.5*Vp, Vp + 0.5*Vp])
    grid on;
    
    title('Moduladora y portadora');
    % legend(["Señal moduladora", "", "", "Señal portadora", "", ""])
    set(gca,'xtick',[]) % Do not plot numbers on x axis
    xlim([0, t(end)])
    grid on;
end

function plot_pwm(t, pwm)
    V_dc_up = max(pwm);
    V_dc_low = min(pwm);

    plot(t,pwm,'k', 'linewidth', 0.75)
    hold on
    plot([0, t(end)], [0, 0], 'k--', 'linewidth', 0.5);
%     set(gca, 'ytick', [V_dc_low, V_dc_up],'yticklabel',...
%         ["V_{-}"; "V_{+}"])

    xlabel('Tiempo');
    title('PWM');
    set(gca,'xtick',[]) % Do not plot numbers on x axis
    ylim([V_dc_low-0.5*V_dc_up, V_dc_up + 0.5*V_dc_up])
    xlim([0, t(end)])
    grid on;

end

function plot_amp_spectrum(hrmncs, min_for_plot)
    P1 = abs(hrmncs);
    numb_hmncs = length(P1);    
    % Horizontal line
    semilogy([min_for_plot, numb_hmncs],...
             [min_for_plot, min_for_plot], 'b', "lineWidth", 1); 
    hold on
    for h = 1:numb_hmncs % Vertical line for each harmonic
        if P1(h) > min_for_plot
            % the index h_index corresponds to the h - 1 harmonic, as the 
            % DC component is at h_index = 1
            semilogy([h-1, h-1], [min_for_plot, P1(h)], 'b', "lineWidth", 1);
        end
    end
    xlabel('Número (f/fm)')
    ylabel('Amplitud')
    if (max(P1) > min_for_plot)
        ylim([0.9*min_for_plot, 1.2*max(P1)])
    end
end




