function sideband_amps = get_region_amps(vm, carrier_type, m)
    % GET_REGION_AMPS Computes sideband amplitudes for a given region.
    % Warning: it does not take into account the effect of neighboring
    % regions, which can be noticeable for large m or small mf.
    %
    % Arguments:
    %   vm (vector): A vector containing one period of the modulating signal.
    %   carrier_type (int): Specifies the type of carrier:
    %                       1 for symmetric carriers,
    %                       2 for sawtooth carriers.
    %   m (int): The region to be computed.
    %
    % Returns:
    %   sideband_amps (vector): Sideband amplitudes.

    if m == 0
        f = 2*vm;
    elseif carrier_type == 1
        f = (4/(m*pi))*sin(m*pi*(1 + vm)/2);
    elseif carrier_type == 2
        f = (1i*2/(m*pi))*(cos(m*pi*vm) - ((-1)^m)- 1i*sin(m*pi*vm));
    else
        error('Invalid carrier_type. It must be 1 (symmetric) or 2 (sawtooth).');
    end
    abs_ifft = abs(ifft(f));
    sideband_amps = abs_ifft(1:round(end/2)); 
end
