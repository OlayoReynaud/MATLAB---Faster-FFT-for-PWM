clear; clc; close all;

% Time vector
t = linspace(0, 2*pi, 5000);

% Sine wave
mod = 0.99*sin(t);

% Triangular wave (carrier)
tri_freq = 30; % Frequency multiplier for triangle wave

numb_levels = 5;
cars = zeros(length(t), numb_levels);
pwms = zeros(length(t), numb_levels);
car_amp =1/numb_levels;

tld_lyot = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
hold on;
multilevel_pwm = zeros(length(t), 1);
for i = 1:numb_levels
    car_dc = car_amp*(-numb_levels-1+2*i);
    cars(:, i) = car_amp*sawtooth(tri_freq * t, 0.5) + car_dc;
    pwms(:, i) = generate_pwm(mod, cars(:, i));
    plot(t, cars(:, i), 'k', 'LineWidth', 1);
end

plot(t, mod, 'k', 'LineWidth', 2);
ylim([-1.1, 1.1]);
xlim([0, 2*pi]);
xlabel('')
ylabel('')
hide_tics()
grid on;
hold off;

nexttile
plot(t, pwms(:, 4), 'k', 'LineWidth', 1);
ylim([-1.1, 1.1]);
xlim([0, 2*pi]);
hide_tics()
xlabel('')
ylabel('')
grid on;


figure;
hold on;
for i = 1:4
    if rem(numb_levels, 2)==0
        plot(mod*(numb_levels-1) + 2*(ceil((numb_levels-1)/2) - i))
    else
        plot(1 + mod*(numb_levels-1) + 2*(ceil((numb_levels-1)/2) - i))
    end
end

function pwm = generate_pwm(mod, car)
    N = length(mod);
    pwm = ones(N, 1);
    for i = 1:N
        if mod(i) < car(i)
            pwm(i) = -1;
        end
    end
end

function hide_tics()
    % Hide tick labels and marks
    ax = gca;
    ax.XColor = 'none';
    ax.YColor = 'none';
    
    % Add grid manually
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    
    % Overlay invisible axes for grid with ticks
    uistack(ax, 'top')
end