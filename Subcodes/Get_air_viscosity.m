function eta_air = Get_air_viscosity(T)
    % https://resources.system-analysis.cadence.com/blog/msa2022-the-relationship-between-the-kinematic-viscosity-of-air-and-temperature
    eta_25_deg = 1.57e-5; % [m^2 s-1]
    eta_air = eta_25_deg / (298^(3/2)) *  T^(3/2);
end