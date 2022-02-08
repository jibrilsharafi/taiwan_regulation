function [E, P_AC] = simulationEdReg(f, Pps_AC, OpMode, E_max, P_nominal, eta, EL_target_0, deadband)
    %% USEFUL VARIABLES
    T = length(f);  % number of seonds (=timesteps)
    dt=1/3600;      % s -> h
    %% OPERATION MODE       
    % Different operation ranges given by TSO
    if strcmp(OpMode,"dReg0.25")
        % dReg 0.25 parameters
        % Map of prescribed response
        upcurve=[100 100 52 9 9 -52 -100 -100];
        lowcurve=[100 100 52 -9 -9 -52 -100 -100];
        freq=[59 59.75 59.86 59.98 60.02 60.14 60.25 61];
    else
        % dReg 0.5 parameters (valid also for EdReg)
        % Map of prescribed response
        % First and last values will be cut when plotting otherwise not right
        % shape
        upcurve=[100 100 48 9 9 -48 -100 -100];
        lowcurve=[100 100 48 -9 -9 -48 -100 -100];
        freq=[59 59.5 59.75 59.98 60.02 60.25 60.5 61];
    end
    
    % Creating the curves that describe the upper and lower bound of the FR.
    % This means that (eg) if it is charging (negative power), the response
    % will be a point (depending on the frequency, that's why interp1) on the
    % lower curve Discharge-oriented response
    dchResponse = interp1(freq, upcurve, f) * P_nominal / 100;
    % Charge-oriented response
    chResponse = interp1(freq, lowcurve, f) * P_nominal / 100;
    
    %% ZERO RESPONSE ARRAY
    % Creating zero response
    % First, we check the positions where NOT both the charging curve and
    % discharging curve (given by operation area upper and lower boundaries)
    % are positive or negative, meaning we only get the middle part with 0%
    % power
    noRespFlag = sign(dchResponse.*chResponse) == -1;
    % We create this array of infinite values as later we will use the min
    % function. The values which will be 0 will be substituted with 0 in the
    % array, the rest stays infinite
    noRespSet = inf(T,1);
    % The values which could actually be 0 power are substituted in
    noRespSet(noRespFlag) = 0;
    % This variable represents the zero power when it is possible, e.g. when
    % the upper bound and the lower bound are one positive and one negative,
    % while it has an inf value when upper and lower bound have the same sign
    % (both above or below the zero power). This ensures that we can choose the
    % best approach (zero power) when we can, while still following the
    % constraints.
    zeroResponse = min(abs([dchResponse'; chResponse'; noRespSet'])) .* sign(60-f');
    % Lowest allowed response (closest to 0% power, to achieve less usage)
    
    %% CREATING EMPTY ARRAYS
    E = nan(1,T+1);
    EL_target = nan(1,T+1);
    EL_up = nan(1,T+1);
    EL_down = nan(1,T+1);
    P_target_AFR_AC = nan(1,T);
    P_AC = nan(1,T);
    P_SOC_managment_AC = zeros(1,T);
    P_SOC_managment_DC = zeros(1,T);
    P_PS_DC = nan(1,T);
    P_DC = nan(1,T);
    P_target_AC = nan(1,T);
    
    % DISPATCH SIMULATION
    % hflag signals if the target power request (energy) is above, below or
    % inside the operational area
    hFlag = nan(1,T+1)';
    hFlag(1) = 0;
    
    %% INITIALISING SOME VALUES
    E(1) = E_max * 0.5;        % Starting SOC: 50%
    EL_target(1) = EL_target_0;
    EL_up(1) = EL_target_0 + deadband;
    EL_down(1) = EL_target_0 - deadband;
    
    %% MAIN LOOP
    for t = 1:T
        % Adjusting reference level for AFR SoC management based on the
        % integration of the peak shaving setpoint (hysteresis follows the
        % curve due to peak shaving)
        EL_up(t) = EL_target(t) + deadband;
        EL_down(t) = EL_target(t) - deadband;
        if EL_up(t) > 0.99
            EL_up(t) = EL_up(t-1);
        end
        if EL_up(t) < 0.01
            EL_down(t) = EL_down(t-1);
        end
    
        % Pmin and Pmax are related to the minimum and maximum constant power
        % available for 1 full second
        % Pmin is the maximum power while charging
        Pmin = -(E_max - E(t)) / (eta*dt);
        % Pmax is the maximum power while discharging
        Pmax = E(t) * (eta/dt);
    
        if hFlag(t) == 0       % Freedom of choice of power
            if Pps_AC(t) == 0       % No peak shaving required -> Do the least possible
                P_target_AFR_AC(t) = zeroResponse(t);
            elseif Pps_AC(t) > 0    % Peak shaving (discharging) -> counteract with opposite (charging)
                P_target_AFR_AC(t) = chResponse(t);
                P_SOC_managment_AC(t) = chResponse(t) - zeroResponse(t);
            elseif Pps_AC(t) < 0    % Peak shaving (charging) -> counteract with opposite (discharging)
                P_target_AFR_AC(t) = dchResponse(t);
                P_SOC_managment_AC(t) = dchResponse(t) - zeroResponse(t);
            end
        elseif hFlag(t) == 1       % Force the upper part of the operational curve
            P_target_AFR_AC(t) = dchResponse(t);
        elseif hFlag(t) == -1      % Force the lower part of the operational curve
            P_target_AFR_AC(t) = chResponse(t);
        end
    
        % Now we calculate the single DC components, then we sum them. This
        % process is needed in order to set the EL_target.
        P_PS_DC(t) = max(Pps_AC(t)/eta, Pps_AC(t)*eta);
        P_SOC_managment_DC(t) = max(P_SOC_managment_AC(t)/eta, P_SOC_managment_AC(t)*eta);
        P_DC(t) = P_PS_DC(t) + P_SOC_managment_DC(t);
    
        % Real power required by grid, AC side
        P_target_AC(t) = P_target_AFR_AC(t) + Pps_AC(t);
        P_AC(t) = min([max([P_target_AC(t), Pmin, -P_nominal]), Pmax, P_nominal]);
    
        % Calculation of energy and energy target for next loop
        E(t+1) = E(t) - max(P_AC(t)/eta, P_AC(t)*eta) * dt;
        EL_target(t+1) = EL_target(t) - P_DC(t) * dt / E_max;
    
        % Depending on the energy level, hFlag decides what to do next loop
        if ((E(t+1)/E_max) > EL_up(t)) && (Pps_AC(t) <= 0)
            hFlag(t) = 1;
        elseif E(t+1)/E_max < EL_down(t) && (Pps_AC(t) >= 0)
            hFlag(t) = -1;
        elseif (((E(t+1)/E_max) < EL_target(t)) && hFlag(t)==1) || (((E(t+1)/E_max) > EL_target(t)) && hFlag(t)==-1)
            hFlag(t) = 0;
        end
        hFlag(t+1) = hFlag(t);
    end
end