clear; clc; close all;

parentDir = pwd; % Parent directory
p1 = fileparts(parentDir);  % One level up
p2 = fileparts(p1);  % One level up
addpath(fullfile(p2, 'auxiliary_functions/'));
addpath(fullfile(p2, 'algorithm_functions/'));

zss_types = ["SPWM", "THIPWM1/6", "CB-SVPWM", "THIPWM1/4", ...
    "SDPWM", "DPWM0"];
numb_zss_types = length(zss_types);


numb_mf = 5e3;
% Get evenly spaced numbers in a logarithmic scale
mf_vals = unique(round(logspace(log10(10), log10(1e5), numb_mf)));
% mf_vals = 10:1:5e4;

numb_mf_vals = length(mf_vals)

colors = lines(numb_zss_types); % Get numb_zss_types distinct colors

alpha = 0.4;
background = [1 1 1];  % white
light_colors = alpha*colors + (1 - alpha) * background; % More transparent


for carrier_type = [1, 2]
    disp_hline(); disp_hline(); disp_hline();
    figure();
    for i_zss = 1:numb_zss_types
        disp_hline(); disp_hline(); disp(zss_types(i_zss))
        pts = zeros(3, numb_mf_vals);
        times = zeros(3, numb_mf_vals);
        errors = zeros(3, numb_mf_vals);
        for i_mf = 1:numb_mf_vals
            mf = mf_vals(i_mf);
            if rem(i_mf, 10) == 0
                disp(strcat(num2str(100*i_mf/numb_mf_vals), " %"))
            end
            [pts(:, i_mf), times(:, i_mf), errors(:, i_mf)] = ...
                test_result(i_zss-1, mf, carrier_type, false);
        end
        % Remove NaNs (somehow it appears once, for mf=15, DPWMO, Sawtooth)
        idx_no_NaNs = ~isnan(pts(3, :));

        % Plot the real values with trasparent markers
        subplot(1, 2, 1); hold on;
        h1(1) = plot(mf_vals(idx_no_NaNs), pts(1, idx_no_NaNs), "*",'Color', ...
            light_colors(i_zss, :));
        h1(2) = plot(mf_vals(idx_no_NaNs), pts(3, idx_no_NaNs), "o", 'Color', ...
            light_colors(i_zss, :));
        set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
        yticks = 2.^(0:24);
        set(gca, 'YTick', yticks);
        set(gca, 'YTickLabel', arrayfun(@(x) ['2^{' num2str(log2(x)) '}'], ...
            yticks, 'UniformOutput', false));
        
        uistack(h1, 'bottom'); % Move trasparent lines to background

        subplot(1, 2, 2); hold on;
        h2(1) = plot(mf_vals(idx_no_NaNs), 1e3*times(1, idx_no_NaNs), "*", ...
            'Color', light_colors(i_zss, :));
        h2(2) = plot(mf_vals(idx_no_NaNs), 1e3*times(3, idx_no_NaNs), "o",...
            'Color', light_colors(i_zss, :));
        set(gca, 'XScale', 'log');
        
        uistack(h2, 'bottom'); % Move trasparent lines to background

        % Plot the mooving averages on top
        k_movmean = 50;
        subplot(1, 2, 1); hold on;
        plot(mf_vals(idx_no_NaNs), ...
            variable_moving_average(mf_vals(idx_no_NaNs), ...
            pts(1, idx_no_NaNs)), "--", 'Color', colors(i_zss, :), ...
            "LineWidth",1.5)
        plot(mf_vals(idx_no_NaNs), ...
            variable_moving_average(mf_vals(idx_no_NaNs), ...
            pts(3, idx_no_NaNs)), "-", 'Color', colors(i_zss, :), ...
            "LineWidth",1.5)

        subplot(1, 2, 2); hold on;
        plot(mf_vals(idx_no_NaNs), ...
            variable_moving_average(mf_vals(idx_no_NaNs), ...
            1e3*times(1, idx_no_NaNs)), "--", 'Color', colors(i_zss, :), ...
            "LineWidth",1.5)
        plot(mf_vals(idx_no_NaNs), ...
            variable_moving_average(mf_vals(idx_no_NaNs), ...
            1e3*times(3, idx_no_NaNs)), "-", 'Color', colors(i_zss, :), ...
            "LineWidth",1.5)

        % % Plot the mooving averages on top
        % k_movmean = 50;
        % subplot(1, 2, 1); hold on;
        % plot(mf_vals(idx_no_NaNs), ...
        %     movmean(pts(1, idx_no_NaNs), k_movmean), "--", 'Color', ...
        %     colors(i_zss, :), "LineWidth",1.5)
        % plot(mf_vals(idx_no_NaNs), ...
        %     movmean(pts(3, idx_no_NaNs), k_movmean), "-", 'Color', ...
        %     colors(i_zss, :), "LineWidth",1.5)
        % 
        % subplot(1, 2, 2); hold on;
        % plot(mf_vals(idx_no_NaNs), ...
        %     movmean(1e3*times(1, idx_no_NaNs), k_movmean), "--", ...
        %     'Color', colors(i_zss, :), "LineWidth",1.5)
        % plot(mf_vals(idx_no_NaNs), ...
        %     movmean(1e3*times(3, idx_no_NaNs), k_movmean), "-",...
        %     'Color', colors(i_zss, :), "LineWidth",1.5)

        
    end
    subplot(1, 2, 1); title("Number of points needed"); 
    xlabel("m_f"); ylabel("Sampling points"); grid on;
    subplot(1, 2, 2); title("Time employed");
    xlabel("m_f"); ylabel("Time (ms)"); grid on;

    for i_subplot = [1, 2] % Create legend
        subplot(1, 2, i_subplot)
        for i_leg1 = 1:numb_zss_types % Dummy plots for color legend
            h_colors(i_leg1) = plot(NaN, NaN, 's', 'Color', ...
                colors(i_leg1, :), 'MarkerFaceColor', colors(i_leg1, :), ...
                'LineWidth', 1.5, 'DisplayName', zss_types(i_leg1));
        end
        h_blank = plot(NaN, NaN, 'w-','LineWidth', 1.5, 'DisplayName', '');
        h_solid = plot(NaN, NaN, 'k-', ...
                'LineWidth', 1.5, 'DisplayName', 'FFT');
        h_dashed = plot(NaN, NaN, 'k--', ...
                'LineWidth', 1.5, 'DisplayName', 'Proposed method');
        leg = legend([h_colors,h_blank, h_solid, h_dashed], 'Location',...
            'northwest');
    end
end

function y_avg = variable_moving_average(x, y_raw)
    % Performs a moving average on y_raw with a variable window size
    % depending on the corresponding value of x.

    N = length(y_raw);
    y_avg = zeros(size(y_raw));

    for i = 1:N
        % half_window = floor(20 * log10(x(i)));
        half_window = 0.1*x(i);

        % Ensure the window is at least 1
        half_window = max(1, half_window);

        % Compute bounds
        lb = max(x(1), x(i) - half_window);
        ub = min(x(end), x(i) + half_window);

        % Find indices of x within bounds
        window_idxs = x >= lb & x <= ub;

        % Compute average
        y_avg(i) = mean(y_raw(window_idxs));
    end
end

function disp_hline()
    disp("------------------------------------------")
end

function [pts_used, time_needed, max_error] = test_result(zss_type, mf, ...
    carrier_type, display_results)
    % Returns the performace of the fft algorithm and our method.
    %   All outputs are is 3x1. The first value coresponds to the new 
    %       algorithm with approximate computation, the second one to the
    %       new algorithm with acute (inlcuding neigboring regions) 
    %       computation and the third one to the fft. 

    % carrier_type is 1 for symmetric carriers and 2 for sawtooth
    max_error = 1e-3;
    numb_sidebands = min(mf, 10);
    m = 1;
    ma = 1;
    % Getting first accurate results with our algorithm and 2^15 sampling
    % points. 
    vect_size_reference = 2^20;
    amps_ref_simple = get_sideband_amps(zss_type, vect_size_reference, ...
        m, numb_sidebands, ma, carrier_type);
    % amps_ref_full includes the effect of neigboring regions.
    amps_ref_full = get_sideband_amps(zss_type, vect_size_reference, m,...
                                     numb_sidebands, ma, carrier_type, mf);

    vect_size_ours1= 256;
    error_ours1_ok = false;
    while ~error_ours1_ok
        tic
        amps_ours1 = get_sideband_amps(zss_type, vect_size_ours1, m,...
                                  numb_sidebands, ma, carrier_type);
        time_ours1 = toc;
        error_ours1 = max(abs(amps_ref_simple - amps_ours1));
        if (error_ours1 > max_error)
            vect_size_ours1 = vect_size_ours1*2;
        else
            error_ours1_ok = true;
        end
    end

    vect_size_ours2= 256;
    error_ours2_ok = false;
    while ~error_ours2_ok
        tic    
        amps_ours2 = get_sideband_amps(zss_type, vect_size_ours2, m,...
                                      numb_sidebands, ma, carrier_type, mf);
        error_ours2 = max(abs(amps_ref_full - amps_ours2));
        time_ours2 = toc;
        if (error_ours2 > max_error)
            vect_size_ours2 = vect_size_ours2*2;
        else
            error_ours2_ok = true;
        end
    end
   
   
    % Initial vectors size
    vect_size_fft = max(2^15, 2^(nextpow2(mf)+1));
    max_vect_size = 2^25;
    error_fft_ok = false;
    while ~error_fft_ok

        vm_fft = generate_vm(ma, zss_type, vect_size_fft);
        tic
        % Time vector normalized to one period of vm
        t_normalized = linspace(0, 1, vect_size_fft);
        if carrier_type == 1
            car = sawtooth(2*pi*mf*t_normalized, 1/2);
        elseif carrier_type == 2
            car = sawtooth(2*pi*mf*t_normalized - pi, 0);
        end
        pwm_sig = generate_pwm(vm_fft, car);
        amps_spectrum_fft = amp_spectrum_fft(pwm_sig);
        % We take the ones on the left because they are less disturbed
        % by the next region.

        amps_fft = flip(amps_spectrum_fft(mf-numb_sidebands+2:mf+1));
        time_fft = toc;
        error_fft = max(abs(amps_ref_full - amps_fft));
        
        if (error_fft > max_error)
            if (vect_size_fft >= max_vect_size)
                % figure()
                % plot(pwm_sig); hold on; plot(vm_fft); plot(car);
                % figure()
                % plot(amps_ref_full); hold on; plot(amps_fft)
                % abs(amps_ref_full - amps_fft)
                vect_size_fft = NaN; time_fft = NaN; error_fft = NaN;
                warning("FFT did not converge to expected amplitudes." + ...
                    "Retutning NaN.")
                break;
            end
            vect_size_fft = vect_size_fft*2;
        else
            error_fft_ok = true;
        end
    end
    
    if display_results
        disp("Sampling points needed our algorithm - approximate computation")
        disp(strcat("2^", num2str(log2(vect_size_ours1))))
        disp("Sampling points needed our algorithm - accurate computation")
        disp(strcat("2^", num2str(log2(vect_size_ours2))))
        disp("Sampling points needed fft")
        disp(strcat("2^", num2str(log2(vect_size_fft))))
        
        disp("Time with our algorithm  - approximate computation")
        disp(strcat(num2str(time_ours1*1e3), " ms"))
        disp("Time with our algorithm - accurate computation")
        disp(strcat(num2str(time_ours2*1e3), " ms")) 
        disp("Time for fft")
        disp(strcat(num2str(time_fft*1e3), " ms"))
        
        disp("Maximum error with our algorithm - approximate computation")
        error_ours1 = max(abs(amps_ref_full - amps_ours1));
        disp("(neighboring regions are used to compute this error)")
        disp(strcat(num2str((1e5)*error_ours1), "e-5"))
        disp("Maximum error with our algorithm - accurate computation")
        disp(strcat(num2str((1e5)*error_ours2), "e-5"))
        disp("Maximum error with fft")
        disp(strcat(num2str((1e5)*error_fft), "e-5"))
    end

    pts_used = [vect_size_ours1; vect_size_ours2; vect_size_fft];
    time_needed = [time_ours1; time_ours2; time_fft];
    max_error = [error_ours1; error_ours2; error_fft];

end

function amps = get_sideband_amps(zss_type, vect_size, m, numb_sidebands,...
                                  ma, carrier_type, mf)
    % Uses our algorithm to get the amplitude of the numb_sidebands first
    % sidebands.
    %     zss_type: type of modulation (0 for SPWM, 1 for THIPWM1/6, 2 for CB-SVPWM,
    %               3 for THIPWM1/4, 4 for SDPWM, 5 for DPWM).
    %     mf (optional): frequency modulation index. If passed, the
    %                   the neigboring sidebands are taken into account
    %                   (they might have a small effect on the amplitudes).
    
    vm = generate_vm(ma, zss_type, vect_size);
    
    if nargin > 6  % Checks if 'mf' is passed as argument
        if mf < 10
            nargin_right = 100;
        elseif mf < 50
            nargin_right = 20;
        elseif mf < 500
            nargin_right = 5;
        else
            nargin_right = 1;
        end
        pwm_hmncs = get_spectrum(vm, mf, carrier_type, m+nargin_right);
        % We take the ones on the left because are less disturbed
        % by the next regions.
        amps = abs(flip(pwm_hmncs(mf-numb_sidebands+2:mf+1)));        
    else
        amps = get_region_amps(vm, carrier_type, m);
        % We are only interested in the first numb_sidebands sidebands
        amps = amps(1:numb_sidebands);
    end
end

function mod_a = generate_vm(ma, zss_type, vect_size)
    %     zss_type: type of modulation (0 for SPWM, 1 for THIPWM1/6, 2 for CB-SVPWM,
    %               3 for THIPWM1/4, 4 for SDPWM, 5 for DPWM).

    %anlge_step = 2*pi/vect_size
    % angle_vect = 0:anlge_step:2 * pi - anlge_step;
    
    angle_vect = linspace(0, 2*pi, vect_size + 1);
    angle_step = angle_vect(2);
    angle_vect(end) = []; % Remove the last element

    angle_vect = angle_vect(1:end-1);

    va_raw = ma * cos(angle_vect - angle_step);

    phase_shift = round((120 / 360) * length(va_raw));
    vb_raw = circshift(va_raw, phase_shift);
    vc_raw = circshift(vb_raw, phase_shift);
    zss = create_zss_from_mods(va_raw, vb_raw, vc_raw, ma, 0, zss_type, 0);
    mod_a = va_raw + zss;    
end

function zss = create_zss_from_mods(va_raw, vb_raw, vc_raw, ma, delta, ...
        zss_type, dpwm_type)
    % Returns the ZSS signal for carrier-based PWM with homopolar injection.
    % Its amplitude is normalized with respect to the amplitude of the carrier.
    %
    % Args:
    %     va_raw, vb_raw, vc_raw: raw (sinusoidal) modulating signals.
    %     ma: amplitude modulation index.
    %     delta: phase of modulation signal.
    %     zss_type: type of modulation (0 for SPWM, 1 for THIPWM1/6, 2 for CB-SVPWM,
    %               3 for THIPWM1/4, 4 for SDPWM, 5 for DPWM).
    %     dpwm_type: type of DPWM, used only if zss_type is DPWM (0 to 6).
    %     psi: angle used in GDPWM (only for DPWM).
    %
    % Returns:
    %     zss: discretized modulating signal over one full period.

    angle_vect = linspace(0, 2*pi, length(va_raw) + 1);
    angle_vect(end) = []; % Remove the last element

    if zss_type == 0  % SPWM
        zss = zeros(1, length(va_raw));        
    elseif zss_type == 1  % THIPWM1/6
        zss = nth_hrmnc_lst(angle_vect, delta, 3, -ma * (1/6));        
    elseif zss_type == 2  % CB-SVPWM
        zss = -0.5 * (max([va_raw; vb_raw; vc_raw]) + min([va_raw; vb_raw; vc_raw]));        
    elseif zss_type == 3  % THIPWM1/4
        zss = nth_hrmnc_lst(angle_vect, delta, 3, -ma * (1/4));        
    elseif zss_type == 4  % SDPWM
        zss = nth_hrmnc_lst(angle_vect, delta, 3, -ma / (2 * pi)) + ...
              nth_hrmnc_lst(angle_vect, delta, 9, -ma / (60 * pi)) + ...
              nth_hrmnc_lst(angle_vect, delta, 15, -ma / (280 * pi));
        
    else  % DPWM
        zss = zzs_for_dpwm(va_raw, vb_raw, vc_raw, dpwm_type);
    end
end

function result = nth_hrmnc_lst(angle_lst, delta, n, amp_nth_hmnc)
    % Computes a discretized cosine wave with frequency 1/(n*2*pi).
    %
    % Args:
    %     angle_lst: A vector of discretized angle values.
    %     delta: Phase shift.
    %     n: Harmonic number. If angle_lst goes from 0 to 2*pi, then n
    %        periods are present in the output.
    %     amp_nth_hmnc: Amplitude of the output wave (default is 1).
    %
    % Returns:
    %     result: A vector of discretized values of a cosine wave with
    %     frequency 2*pi*n.
    
    if nargin < 4
        amp_nth_hmnc = 1; % Default amplitude
    end
    
    result = amp_nth_hmnc * cos(n * (angle_lst - delta));
end

function zss = zss_function_dpwm(va_raw, vb_raw, vc_raw, condsa, condsb)
    % Compute the homopolar (zss) signal for DPWM.
    % This function evaluates the conditions to determine the phase used and computes:
    % zss = (sign(vi) - vi) based on the conditions.
    %
    % Args:
    %     va_raw, vb_raw, vc_raw: raw modulating signals (before adding the zss).
    %     condsa, condsb: conditions for phase selection (logical arrays).
    %
    % Returns:
    %     zss: discretized zero-sequence-signal

    lists_length = length(va_raw);
    zss = zeros(1, lists_length);  % Initialization

    for i = 1:lists_length
        if condsa(i)
            zss(i) = sign(va_raw(i)) - va_raw(i);
        elseif condsb(i)
            zss(i) = sign(vb_raw(i)) - vb_raw(i);
        else
            zss(i) = sign(vc_raw(i)) - vc_raw(i);
        end
    end
end

function zss = zzs_for_dpwm(va_raw, vb_raw, vc_raw, dpwm_type, psi)
    % Compute the homopolar (zss) signal for different types of DPWM.
    %
    % Args:
    %     va_raw, vb_raw, vc_raw: raw modulating signals.
    %     dpwm_type: type of DPWM.
    %     psi: angle used in GDPWM (default 30 degrees).
    %
    % Returns:
    %     zss: discretized zero-sequence-signal

    if nargin < 5
        psi = 30.0;
    end

    VECT_LENGTH = length(va_raw);

    if dpwm_type == 5  % DPWMMax
        condsa = (va_raw > vb_raw) & (va_raw > vc_raw);
        condsb = (vb_raw > va_raw) & (vb_raw > vc_raw);
        zss = zss_function_dpwm(va_raw, vb_raw, vc_raw, condsa, condsb);

    elseif dpwm_type == 6  % DPWMMin
        condsa = (va_raw < vb_raw) & (va_raw < vc_raw);
        condsb = (vb_raw < va_raw) & (vb_raw < vc_raw);
        zss = zss_function_dpwm(va_raw, vb_raw, vc_raw, condsa, condsb);

    else
        if dpwm_type == 0      % DPWM0
            roll_val = round((330/360) * VECT_LENGTH);
        elseif dpwm_type == 2  % DPWM2
            roll_val = round((30/360) * VECT_LENGTH);
        elseif dpwm_type == 4  % GDPWM
            roll_val = round(((psi - 30) / 360) * VECT_LENGTH);
        else                   % DPWM1 and DPWM3
            roll_val = 0;
        end

        vax = circshift(va_raw, roll_val);
        vbx = circshift(vb_raw, roll_val);
        vcx = circshift(vc_raw, roll_val);

        abs_vax = abs(vax);
        abs_vbx = abs(vbx);
        abs_vcx = abs(vcx);

        if dpwm_type == 3  % DPWM3
            condsa = ((abs_vbx > abs_vax) & (abs_vax > abs_vcx)) | ...
                     ((abs_vcx > abs_vax) & (abs_vax > abs_vbx));
            condsb = ((abs_vax > abs_vbx) & (abs_vbx > abs_vcx)) | ...
                     ((abs_vcx > abs_vbx) & (abs_vbx > abs_vax));
        else
            condsa = (abs_vax > abs_vbx) & (abs_vax > abs_vcx);
            condsb = (abs_vbx > abs_vax) & (abs_vbx > abs_vcx);
        end

        zss = zss_function_dpwm(va_raw, vb_raw, vc_raw, condsa, condsb);
    end
end


