classdef Calc_Gaussian_Rayleigh < handle
    % Useful expressions https://www.mathworks.com/help/matlab/ref/mustbepositive.html
    % https://www.mathworks.com/help/matlab/ref/arguments.html
    properties (SetAccess = public, Hidden = false)
%         NA_det {mustBePositive} = 0.8;  
%         use_GPU logical = true;             %GPU acceleration

        utility;
        constants;

        pGauss; % Gaussian optics parameters
        pOptomech; % Optomechanical parameters
    end

    methods
    % Main function

        function h=Calc_Gaussian_Rayleigh(util) 

        % Physical constants
            h.constants = struct;
            h.constants.epsilon0 =  8.8541878128e-12; % [F⋅m−1] Vacuum permittivity
            h.constants.c = 299792458; % [m s-1] Light speed
            h.constants.m_gas = 28.964/(6.02214076*10^23)*10^-3; % [kg; N2 gas mass]
            h.constants.kB = 1.380649e-23; % [m^2 kg s-2 K-1] Boltzmann constant

        % Initialize utility
            utility = struct;
            % Optics parameters
            utility.lambda = 1064e-9; % [m]
            utility.NA_foc = 0.8;
            utility.NA_det = 0.8;
            utility.Power = 100e-3; 
            % Sample parameters
            utility.radius = 69e-9; % [m]
            utility.density = 2.196e3;
            utility.RI = GET_RI_SiO2(utility.lambda);
            % Environmental damping parameter
            utility.p = 6.3e2; % [Pa; 1 Bar = 10^5 Pa; 6.3 mBar = 630 Pa]
            utility.T = 298; % [K; temperature]
            % GPU parameters
            utility.use_GPU = false;    

        % Update parameters
            h.utility = utility;
            if nargin >= 1
                h.utility=SET_update_params(h.utility,util);
            end

        % Set physical parameteres according to the input values
        % Particle volume & mass
            h.utility.k = 2*pi / h.utility.lambda;
            h.utility.Vol = 4*pi*(h.utility.radius)^3 / 3; % [m^3]
            h.utility.Mass = h.utility.Vol * h.utility.density; % [kg]

        % Power-normalized Rayleigh scattering parameters
        % Jan Gieseler PhD dissertation "Dynamics of optically levitated nanoparticles in high vacuum"
        % Chapter 2
            h.utility.alpha = 3 * h.utility.Vol * h.constants.epsilon0 *...
                (h.utility.RI^2-1) / (h.utility.RI^2 + 2); % polarizability [F m^2]
            h.utility.alpha_eff = h.utility.alpha / ... % Effective polarizability due to radiation reaction, Eq (2.13)
                (1 - 1i*(h.utility.k^3/6/pi/h.constants.epsilon0)*h.utility.alpha);
            h.utility.Sigma_scat = h.utility.k^4 * abs(h.utility.alpha)^2 ./...
            (6*pi*h.constants.epsilon0^2); % Scattering cross section

        % Gas damping parameter
            h.utility.eta_air = Get_air_viscosity(h.utility.T);
            h.utility.l_free = h.utility.eta_air / h.utility.p * sqrt(pi*h.constants.kB*h.utility.T/...
                2/h.constants.m_gas); % Mean free path
            h.utility.Kn = h.utility.l_free / h.utility.radius; % Knudsen number
            h.utility.G_gas = 6*pi*h.utility.eta_air*h.utility.radius / h.utility.Mass *...
                0.619/(0.619+h.utility.Kn)*(1+(0.31*h.utility.Kn/(0.785+1.152*h.utility.Kn+h.utility.Kn^2)));

        end

        function pGauss = Calc_pGauss(h, varargin)
            % Initializer params and parser: https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
            pGauss = struct;
            p = inputParser;

            defaultMode = 'paraxial';
            validMode = {'paraxial', 'empirical', 'vectoral'};
            checkMode = @(x) any(validatestring(x, validMode));
            addOptional(p, 'Mode', defaultMode, checkMode)
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            
        % Get beam waists & z0 (Rayleigh range)
            switch p.Results.Mode 
                case 'paraxial'
                % Mode 1. Paraxial Gaussian optics
                % https://www.rp-photonics.com/numerical_aperture.html
                % https://en.wikipedia.org/wiki/Gaussian_beam
                    pGauss.waistX =   h.utility.lambda / pi / h.utility.NA_foc;
                    pGauss.waistY = pGauss.waistX;
                    pGauss.z0 = pGauss.waistX /  h.utility.NA_foc;
                case 'empirical'
                % Mode 2. Empirical value from Jan Gieseler PhD dissertation
                % "Dynamics of optically levitated nanoparticles in high vacuum"
                % Chapter. 2
                    pGauss.waistX = 687e-9;
                    pGauss.waistY = 542e-9;
                    pGauss.z0 = 1362e-9;
                case 'vectoral' % To be added
            end

            h.pGauss = pGauss;
            h.utility.Intensity = h.utility.Power * 4 / ...
                (pi*h.constants.c*h.constants.epsilon0*...
                pGauss.waistX * pGauss.waistY);
        end

        function pOptomech = Calc_pOptomech(h)

            pOptomech = struct;
        % Trap stiffness
            pOptomech.kX = real(h.utility.alpha_eff) * h.utility.Intensity / h.pGauss.waistX^2;
            pOptomech.kY = real(h.utility.alpha_eff) * h.utility.Intensity / h.pGauss.waistY^2;
            pOptomech.kZ = real(h.utility.alpha_eff) * h.utility.Intensity / 2/h.pGauss.z0^2;

        % Mechanical resonance frequency
            pOptomech.Wx = sqrt(pOptomech.kX/h.utility.Mass);
            pOptomech.Wy = sqrt(pOptomech.kY/h.utility.Mass);
            pOptomech.Wz = sqrt(pOptomech.kZ/h.utility.Mass);


        % Recoil rate - Eqs (3.48) ~ (3.49)
            pOptomech.GbaX = (1/5) *h.utility.Intensity*h.utility.Sigma_scat / h.utility.Mass / ...
                h.constants.c *(2*pi/h.utility.lambda) / pOptomech.Wx * (h.constants.epsilon0*h.constants.c/2);
            pOptomech.GbaY = (2/5) *h.utility.Intensity*h.utility.Sigma_scat / h.utility.Mass / ...
                h.constants.c *(2*pi/h.utility.lambda) / pOptomech.Wy * (h.constants.epsilon0*h.constants.c/2);
            pOptomech.GbaZ = (2/5) *h.utility.Intensity*h.utility.Sigma_scat / h.utility.Mass / ...
                h.constants.c *(2*pi/h.utility.lambda) / pOptomech.Wz * (h.constants.epsilon0*h.constants.c/2);
            h.pOptomech = pOptomech;
        end



    end


end