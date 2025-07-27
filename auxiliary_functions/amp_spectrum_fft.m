function P1 = amp_spectrum_fft(signal)
    % Returns the amplitudes and frequencies of the harmonics using the 
    % fft method. 
    vect_size = length(signal);    
    fft_signal = fft(signal);    
    P2 = abs(fft_signal/vect_size); % Two-sided spectrum
    P1 = P2(1:round(vect_size/2+1)); % One-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);
end
