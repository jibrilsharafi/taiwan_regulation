%% CLEAR EVERYTHING
close all; clear; clc

%% LOAD DATA
% Profiles of frequency and peak shaving setpoint
% Importing 2 variables: f (frequency, mean value 60Hz) and Pps_AC 
% (power peak shaving AC, which has 3 values: -50, 0, 50)
load('ProfiliHoping_NoSoCRef')

% Array of datetime relative to the loaded data. Timestep: seconds
Tst=datetime(2018,12,01,0,0,0):seconds(1):datetime(2018,12,31,23,59,59);
%% STUDY DISTRIBUTION f
% Plot of the residuals and distribution of the frequency
%studyDistribution_f(f);

%% PARAMETERS
% Storage and Service parameters
E_max = 286;            % [MWh]
P_nominal = 100;        % [MW]
eta=0.94;               % [-] - Fixed efficiency

% Initialize values
EL_target_0=0.1;   % Starting target of SOC - 0.1 as we start charging +40%

% Operational mode
OpMode = "EdReg";      %EdReg, dReg05, dReg025

%% INITIALIZING VARIABLES
steps = 50;  % Only a square number so that if fits the subplots evenly
deadband = linspace(0.01, 0.30, steps);

seconds_P_overdch = nan(1, steps);
seconds_P_overch = nan(1, steps);
seconds_E_full = nan(1, steps);
seconds_E_empty = nan(1, steps);

revenue = nan(1, steps);
rev_MW = nan(1, steps);
cycl_num = nan(1, steps);
SoCref_cycl_num = nan(1, steps);

dt=1/3600;      % s -> h

T = length(f);  % Number of seconds

%% MAIN LOOP
f1 = figure("Name", "Frequency response");
for ii = 1:steps
    [E, P_AC] = simulationEdReg(f, Pps_AC, OpMode, E_max, P_nominal, eta, EL_target_0, deadband(ii));

    P_DC = max(P_AC/eta, P_AC*eta);

    % Economic data
    c_el = 1.8;                         % [NTD/kWh]
    dh_capacity_price = 550;            % [NTD/MW/h]
    performance_price = 350;            % [NTD/MW/h]
    net_exchange = sum(P_AC)*dt;           % [MWh]

    revenue(ii) = P_nominal*T*dt*(dh_capacity_price + performance_price) + min(net_exchange,0) * c_el * 1000; % [NTD]

    rev_MW(ii) = revenue(ii) * 365 / (T*dt/24) / P_nominal * 0.03;           % [€]
    cycl_num(ii) = sum(abs(P_DC*dt)) / E_max * 365 / (T*dt/24) / 2;    % [-]

    seconds_E_full(ii) = length(find(E/E_max >= 1));
    seconds_E_empty(ii) = length(find(E/E_max <= 0));
end

%% EVALUATION PLOT
figure()
title('Evaluation performances')
xlabel('Deadband')
legend
hold on

yyaxis left
plot(deadband, cycl_num, 'k', DisplayName='Battery cycles (left axis)', LineWidth=1)

yyaxis right
plot(deadband, seconds_E_full, '-r', DisplayName='Seconds E full', LineWidth=1)
plot(deadband, seconds_E_empty, '-b', DisplayName='Seconds E empty', LineWidth=1)

%% OLD PLOTS
% % Main plot: power&energy
% figure(ii+1)
% hold on
% plot(Tst,P_AC)
% yyaxis right
% plot(Tst,E(1:end-1)/E_max*100, 'r',LineWidth=1)
% hold on
% plot(Tst,EL_target(1:end-1)*100,'--k')
% plot(Tst,EL_up(1:end-1)*100,'--g')
% plot(Tst,EL_down(1:end-1)*100,'--y')
% legend('% power [-]','% energy [-]')
% title(sprintf('Deadband = %0.2f, cycles = %0.2f', [deadband(ii), cycl_num(ii)]))
% hold off
%
% %%
% % Plot sorted power
% figure(2)
% plot(T*dt-(1:T)*dt, (sort(abs(E_cycled_AC))))
% ylabel('Battery power [MW]')
% xlabel('# of operating hours at higher battery power [h]')

%     % Subplot with operational area
%     figure(f1)
%     subplot(sqrt(steps),sqrt(steps),ii)
%     hold on
%     xlabel('Frequency [Hz]')
%     ylabel('% nominal power [-]')
%     title(sprintf('Deadband = %0.2f, cycles = %0.2f', [deadband(ii), cycl_num(ii)]))
%     axis([59.5 60.5 -120 120])
% 
%     OpMode = "EdReg";      %EdReg, dReg05, dReg025
%     
%     % Different operation ranges given by TSO
%     if strcmp(OpMode,"dReg0.25")
%         % dReg 0.25 parameters
%         % Map of prescribed response
%         upcurve=[100 100 52 9 9 -52 -100 -100];
%         lowcurve=[100 100 52 -9 -9 -52 -100 -100];
%         freq=[59 59.75 59.86 59.98 60.02 60.14 60.25 61];
%     else
%         % dReg 0.5 parameters (valid also for EdReg)
%         % Map of prescribed response
%         % First and last values will be cut when plotting otherwise not right
%         % shape
%         upcurve=[100 100 48 9 9 -48 -100 -100];
%         lowcurve=[100 100 48 -9 -9 -48 -100 -100];
%         freq=[59 59.5 59.75 59.98 60.02 60.25 60.5 61];
%     end
%     
%     plot(freq(2:end-1),upcurve(2:end-1),'--k','linewidth',2)
%     plot(freq(2:end-1),lowcurve(2:end-1),'--k','linewidth',2)
% 
%     plot(freq(2:end-1),lowcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
%     plot(freq(2:end-1),upcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
% 
%     plot(freq(2:end-1),lowcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
%     plot(freq(2:end-1),upcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
% 
%     scatter(f(1:10:end), P_AC(1:10:end), 5, E(1:10:end-1)/E_max*100)
%     c = colorbar;
%     c.Label.String = '% energy level [-]';