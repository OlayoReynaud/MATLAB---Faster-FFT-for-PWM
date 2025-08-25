% Get the full path of the current script
currentScriptPath = fileparts(mfilename('fullpath'));
% Get the parent directory
parentDir = fileparts(currentScriptPath);
addpath(parentDir);


K = 1e3;
x = linspace(0,2 * pi * (K - 1)/K,K);
m_f = 40;

T_regular = zeros(K,1);
T_general = zeros(K,7);

for i = 1 : 5e2
    v_mod = zeros(1,K);
    num_harmonics = round(15 * rand);
    for j = 1 : num_harmonics
        v_mod = v_mod + rand/num_harmonics * sin(j * x + rand * 2 * pi);
    end
    % Regular method:
    tic
    fft_for_pwm(v_mod,1,-1,m_f,1,5);
    T_regular(i) = toc;
    for num_levels = 2 : 8   
        % General method:
        carrier_phase = zeros(1,num_levels - 1);
        tic
        pwmfft(v_mod,1,-1,m_f,5,num_levels,carrier_phase);
        T_general(i,num_levels - 1) = toc;
    end

end

disp("Regular method execution time: " + num2str(mean(T_regular)) + " s.");

for i = 1 : 7
    disp("General method execution time (" + num2str(i + 1) + " levels): " + num2str(mean(T_general(:,i))) + " s.");
end
