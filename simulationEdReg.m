function [E, P_AC] = simulationEdReg(deadband, E_max, EL_target_0)
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
    
    E(1) = E_max * 0.5;        % Starting SOC: 50%
    EL_target(1) = EL_target_0;
    EL_up(1) = EL_target_0 + deadband(ii);
    EL_down(1) = EL_target_0 - deadband(ii);
    
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