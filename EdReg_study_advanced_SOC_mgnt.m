%% CLEAR EVERYTHING
close all; clear; clc

%% LOAD DATA

% Profiles including the SoC reference according to "standard" SoC
% management strategy (independent composition of peak shaving & dReg)
%load('ProfiliHoping')

% Profiles of frequency and peak shaving setpoint Importing 2 variables: f
% (frequency, mean value 60Hz) and Pps_AC (power peak shaving AC, which has
% 3 values: -50, 0, 50)
load('ProfiliHoping_NoSoCRef')
%% PARAMETERS
% Array of datetime, timestep: seconds
Tst=datetime(2018,12,01,0,0,0):seconds(1):datetime(2018,12,31,23,59,59);

% Storage and Service parameters
E_max = 286;         % [MWh]
P_nominal = 100;        % [MW]

% Initialize values
EL_target_0=0.1;   % Starting target of SOC - 0.1 as we start charging +40%
eta=0.94;     % Fixed efficiency
dt=1/3600;      % s -> h

%% OPTIONS
deadband=0.05;  %% Hysteresis band

% Test with different datasets 
%f(:) = 60;   % Fixed 60Hz
%f(:) = 60 + 0.05*randn(length(f),1); % 60Hz mean, 0.05 std

%% OPERATION MODE
OpMode = "EdReg";      %EdReg, dReg05, dReg025

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

%% INITIALIZING VARIABLES
% Number of seconds
T = length(dchResponse);

% Empty arrays
E = nan(1,T+1);
EL_target = nan(1,T+1);
EL_up = nan(1,T+1);
EL_down = nan(1,T+1);
P_target_AFR_AC = nan(1,T);
P_AC = nan(1,T);
P_SOC_managment_AC = zeros(1,T);
P_SOC_managment_DC = zeros(1,T);

P_PS_DC = nan(1,T);
P_SOC_PS_AC = nan(1,T);
P_DC = nan(1,T);
P_target_AC = nan(1,T);

% Initial values
E(1) = E_max * 0.5;        % Starting SOC: 50%
EL_target(1) = EL_target_0;
EL_up(1) = EL_target(1) + deadband;
EL_down(1) = EL_target(1) - deadband;

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

%% DISPATCH SIMULATION
% hflag signals if the target power request (energy) is above, below or
% inside the operational area
hFlag = nan(1,T+1)';
hFlag(1) = 0;

for t = 1:T
    % Adjusting reference level for AFR SoC management based on the
    % integration of the peak shaving setpoint (hysteresis follows the
    % curve due to peak shaving)
    EL_up(t) = EL_target(t) + deadband;
    EL_down(t) = EL_target(t) - deadband;

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
    if (E(t+1)/E_max) > EL_up(t)
        hFlag(t) = 1;
    elseif E(t+1)/E_max < EL_down(t)
        hFlag(t) = -1;
    elseif (((E(t+1)/E_max) < EL_target(t)) && hFlag(t)==1) || (((E(t+1)/E_max) > EL_target(t)) && hFlag(t)==-1)
        hFlag(t) = 0;
    end
    hFlag(t+1) = hFlag(t); 
end

E_cycled_AC = max(P_AC/eta, P_AC*eta);

%% ECONOMIC EVALUATION
% Economic data
c_el = 1.8;                         % [NTD/kWh]
dh_capacity_price = 550;            % [NTD/MW/h]
performance_price = 350;            % [NTD/MW/h]
net_exchange = sum(P_AC)*dt;           % [MWh]

revenue = P_nominal*T*dt*(dh_capacity_price + performance_price) + min(net_exchange,0) * c_el * 1000; % [NTD]

rev_MW = revenue * 365 / (T*dt/24) / P_nominal * 0.03;           % [€]
cycl_num = sum(abs(E_cycled_AC*dt)) / E_max * 365 / (T*dt/24) / 2;    % [-]
fprintf('Battery cycles: %f\n', cycl_num)

%% PLOTS
% 2 main subplots: power&energy, frequency
figure(1)
hold on
plot(Tst,P_AC)
yyaxis right 
plot(Tst,E(1:end-1)/E_max*100, 'r',LineWidth=1)
hold on
plot(Tst,EL_target(1:end-1)*100,'--k')
plot(Tst,EL_up(1:end-1)*100,'--g')
plot(Tst,EL_down(1:end-1)*100,'--y')
legend('% power [-]','% energy [-]')
title('Dispatch profile')

% Plot sorted power
figure(2)
plot(T*dt-(1:T)*dt, (sort(abs(E_cycled_AC))))
ylabel('Battery power [MW]')
xlabel('# of operating hours at higher battery power [h]')

% Plot frequency vs power response
figure(3)
hold on
plot(freq(2:end-1),upcurve(2:end-1),'--k','linewidth',2)
plot(freq(2:end-1),lowcurve(2:end-1),'--k','linewidth',2)

plot(freq(2:end-1),lowcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
plot(freq(2:end-1),upcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)

plot(freq(2:end-1),lowcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
plot(freq(2:end-1),upcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)

scatter(f(1:100:end),P_AC(1:100:end), 5, E(1:100:end-1)/E_max,'filled')
colorbar
xlabel('Frequency [Hz]')
ylabel('% nominal power [-]')