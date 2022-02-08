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
% To plot the residuals, we reduce the number of points by a reduction
% factor
reduction_factor = 1e3;
f_reduced = f(1:reduction_factor:end);
% To plot an historam of the distribution of the frequency, we choose the
% number of bins according to Sturge’s Rule
n_bins = floor(1 + 3.322*log10(T));

% Create a subplot
figure()
subplot(1,2,1)
hold on
grid on
xlabel('Frequency [Hz]')
ylabel('Count [-]')
legend(Location="northwest")
histogram(f,NumBins=n_bins,DisplayName=sprintf('Mean = %2.4f, stf = %0.4f', [mean(f), std(f)]))

subplot(1,2,2)
hold on
grid on
xlabel('Time [s]')
ylabel('Residuals [Hz]')
legend
plot(1:reduction_factor:T, f_reduced-mean(f_reduced), 'o', DisplayName='Residuals')
plot([0 T], [0 0], 'r--', LineWidth=1, DisplayName='Mean value')
%% PARAMETERS
% Storage and Service parameters
E_max = 286;         % [MWh]
P_nominal = 100;        % [MW]

% Initialize values
EL_target_0=0.1;   % Starting target of SOC - 0.1 as we start charging +40%
eta=0.94;     % Fixed efficiency
dt=1/3600;      % s -> h

% Number of seconds
T = length(Tst);

%% OPTIONS
%deadband=0.05;  %% Hysteresis band

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

%% INITIALIZING VARIABLES
steps = 4;  % Only a square number so that if fits the subplots evenly
deadband = linspace(0.01, 0.13, steps);

seconds_P_overdch = nan(1, steps);
seconds_P_overch = nan(1, steps);
seconds_E_full = nan(1, steps);
seconds_E_empty = nan(1, steps);

revenue = nan(1, steps);
rev_MW = nan(1, steps);
cycl_num = nan(1, steps);
SoCref_cycl_num = nan(1, steps);

%% MAIN LOOP
for ii = 1:steps

    P_DC_from_AC = max(P_AC/eta, P_AC*eta);

    % Economic data
    c_el = 1.8;                         % [NTD/kWh]
    dh_capacity_price = 550;            % [NTD/MW/h]
    performance_price = 350;            % [NTD/MW/h]
    net_exchange = sum(P_AC)*dt;           % [MWh]

    revenue(ii) = P_nominal*T*dt*(dh_capacity_price + performance_price) + min(net_exchange,0) * c_el * 1000; % [NTD]

    rev_MW(ii) = revenue(ii) * 365 / (T*dt/24) / P_nominal * 0.03;           % [€]
    cycl_num(ii) = sum(abs(P_DC_from_AC*dt)) / E_max * 365 / (T*dt/24) / 2;    % [-]

    index_P_overdch = find(P_AC/P_nominal >= 1);
    index_P_overch = find(P_AC/P_nominal <= -1);
    index_E_full = find(E/E_max >= 1);
    index_E_empty = find(E/E_max <= 0);

    seconds_P_overdch(ii) = length(index_P_overdch);
    seconds_P_overch(ii) = length(index_P_overch);
    seconds_E_full(ii) = length(index_E_full);
    seconds_E_empty(ii) = length(index_E_empty);

    % Subplot with operational area
    figure(1)
    subplot(sqrt(steps),sqrt(steps),ii)
    hold on
    xlabel('Frequency [Hz]')
    ylabel('% nominal power [-]')
    title(sprintf('Deadband = %0.2f, cycles = %0.2f', [deadband(ii), cycl_num(ii)]))
    axis([59.5 60.5 -120 120])

    plot(freq(2:end-1),upcurve(2:end-1),'--k','linewidth',2)
    plot(freq(2:end-1),lowcurve(2:end-1),'--k','linewidth',2)

    plot(freq(2:end-1),lowcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
    plot(freq(2:end-1),upcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)

    plot(freq(2:end-1),lowcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
    plot(freq(2:end-1),upcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)

    scatter(f(1:10:end), P_AC(1:10:end), 5, E(1:10:end-1)/E_max*100)
    c = colorbar;
    c.Label.String = '% energy level [-]';

    % Main plot: power&energy
    figure(ii+1)
    hold on
    plot(Tst,P_AC)
    yyaxis right
    plot(Tst,E(1:end-1)/E_max*100, 'r',LineWidth=1)
    hold on
    plot(Tst,EL_target(1:end-1)*100,'--k')
    plot(Tst,EL_up(1:end-1)*100,'--g')
    plot(Tst,EL_down(1:end-1)*100,'--y')
    legend('% power [-]','% energy [-]')
    title(sprintf('Deadband = %0.2f, cycles = %0.2f', [deadband(ii), cycl_num(ii)]))
    hold off
end

%% EVALUATION PLOT
figure()
title('Evaluation performances')
xlabel('Deadband')
legend
hold on

yyaxis left
plot(deadband, cycl_num, 'k', DisplayName='Battery cycles', LineWidth=1)

yyaxis right
plot(deadband, seconds_E_full, '-r', DisplayName='Seconds E full', LineWidth=1)
plot(deadband, seconds_E_empty, '-b', DisplayName='Seconds E empty', LineWidth=1)

%% OLD PLOTS
% %% PLOTS
% % 2 main subplots: power&energy, frequency
% figure(1)
% hold on
% plot(Tst,P_AC)
% yyaxis right
% plot(Tst,E(1:end-1)/E_max*100, 'r',LineWidth=1)
% hold on
% plot(Tst,EL_target(1:end-1)*100,'--k')
% plot(Tst,EL_up(1:end-1)*100,'--g')
% plot(Tst,EL_down(1:end-1)*100,'--y')
% legend('% power [-]','% energy [-]')
% title('Dispatch profile')
%
% %%
% % Plot sorted power
% figure(2)
% plot(T*dt-(1:T)*dt, (sort(abs(E_cycled_AC))))
% ylabel('Battery power [MW]')
% xlabel('# of operating hours at higher battery power [h]')
%
% %%
% % Plot frequency vs power response
% figure(3)
% hold on
% plot(freq(2:end-1),upcurve(2:end-1),'--k','linewidth',2)
% plot(freq(2:end-1),lowcurve(2:end-1),'--k','linewidth',2)
%
% plot(freq(2:end-1),lowcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
% plot(freq(2:end-1),upcurve(2:end-1)+max(Pps_AC),'--k','linewidth',2)
%
% plot(freq(2:end-1),lowcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
% plot(freq(2:end-1),upcurve(2:end-1)+min(Pps_AC),'--k','linewidth',2)
%
% scatter(f(1:10:end),P_AC(1:10:end), 5, E(1:10:end-1)/E_max*100)
%
% c = colorbar;
% c.Label.String = '% energy level [-]';
% xlabel('Frequency [Hz]')
% ylabel('% nominal power [-]')