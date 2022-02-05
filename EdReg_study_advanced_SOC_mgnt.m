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
Emax = 286;         % [MWh]
Pserv = 100;        % [MW]

% Initialize values
ELtarget=0.1;   % Starting target of SOC
etaOW=0.94;     % Fixed efficiency
dt=1/3600;      % s -> h

%% OPTIONS
deadband=0.05;  %% Hysteresis band

% Test with different datasets 
%f(:) = 60;   % Fixed 60Hz
%f(:) = 60 + 0.03 .*randn(length(f),1); % 60Hz mean, 0.05 std

% Check the width of the peak shaving square wave
% index = find(Pps_AC);
% width = [];
% for i = 1:length(index)-1
%     if (index(i+1) - index(i)) > 1
%         width = [width, (index(i+1) - index(i))*(1+-2*((index(i+1) - index(i)) < 20000))];
%     end
% end

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
dchresponse = interp1(freq, upcurve, f) * Pserv / 100;
% Charge-oriented response
chresponse = interp1(freq, lowcurve, f) * Pserv / 100;

%% INITIALIZING VARIABLES
% Number of seconds
T = length(dchresponse);

% Empty arrays
E = nan(1,T+1);
ELref = nan(1,T+1);
ELup = nan(1,T+1);
ELdwn = nan(1,T+1);
P_pcc_target_AFR = nan(1,T);
P_pcc = nan(1,T);
P_soc_mgnt = zeros(1,T);

Pps_AC_adj = nan(1,T);
Pps_DC_adj = nan(1,T);
P_pcc_target = nan(1,T);

% Initial values
E(1) = Emax * 0.5;        % Starting SOC: 50%
ELref(1) = ELtarget;
ELup(1) = ELtarget + deadband;
ELdwn(1) = ELtarget - deadband;

% Creating zero response
norespflag = sign(dchresponse.*chresponse) == -1;
percnoresp = sum(norespflag) / T;
norespset = inf(T,1);
norespset(norespflag) = 0;   
% This variable represents the zero power when it is possible, e.g. when
% the upper bound and the lower bound are one positive and one negative,
% while it has an inf value when upper and lower bound have the same sign
% (both above or below the zero power). This ensures that we can choose the
% best approach (zero power) when we can, while still following the
% constraints.
zeroresponse = min(abs([dchresponse'; chresponse'; norespset'])) .* sign(60-f');
% Lowest allowed response (closest to 0% power, to achieve less usage)

%% DISPATCH SIMULATION
% hflag signals if the target power request (energy) is above, below or
% inside the operational area
hflag = 0;

for t = 1:T
    if strcmp(OpMode, "EdReg")
        % Adjusting reference level for AFR SoC management based on the
        % integration of the peak shaving setpoint (hysteresis follows the
        % curve due to peak shaving)
        ELup(t) = ELref(t) + deadband;
        ELdwn(t) = ELref(t) - deadband;
    end

    % Pmin and Pmax are related to the minimum and maximum constant power
    % available for 1 full second
    % Pmin is the maximum power while charging
    Pmin = -(Emax - E(t)) / (etaOW*dt);
    % Pmax is the maximum power while discharging
    Pmax = E(t) * (etaOW/dt);
    
    if hflag == 0       % Freedom of choice of power 
        if Pps_AC(t) == 0       % No peak shaving required
            P_pcc_target_AFR(t) = zeroresponse(t);
        elseif Pps_AC(t) > 0    % Peak shaving -> counteract with opposite
            P_pcc_target_AFR(t) = chresponse(t);  
            P_soc_mgnt(t) = chresponse(t) - zeroresponse(t);
        elseif Pps_AC(t) < 0    % Peak shaving -> counteract with opposite
            P_pcc_target_AFR(t) = dchresponse(t);
            P_soc_mgnt(t) = dchresponse(t) - zeroresponse(t);
        end            
    elseif hflag == 1       % Force the upper part of the operational curve
        P_pcc_target_AFR(t) = dchresponse(t);
    elseif hflag == -1      % Force the lower part of the operational curve
        P_pcc_target_AFR(t) = chresponse(t);
    end

    % Power AC side including the SOC managment
    Pps_AC_adj(t) = Pps_AC(t) + P_soc_mgnt(t);
    % Real power DC side 
    Pps_DC_adj(t) = max(Pps_AC_adj(t)/etaOW, Pps_AC_adj(t)*etaOW);
    
    P_pcc_target(t) = P_pcc_target_AFR(t) + Pps_AC(t);
    
    P_pcc(t) = min([max([P_pcc_target(t), Pmin, -Pserv]), Pmax, Pserv]);
    
    E(t+1) = E(t) - max(P_pcc(t)/etaOW, P_pcc(t)*etaOW) * dt;
    ELref(t+1) = ELref(t) - Pps_DC_adj(t) * dt / Emax;

    if (E(t+1)/Emax) > ELup(t)
        hflag = 1;
    elseif E(t+1)/Emax < ELdwn(t)
        hflag = -1;
    elseif (((E(t+1)/Emax) < ELref(t)) && hflag==1) || (((E(t+1)/Emax) > ELref(t)) && hflag==-1)
        hflag = 0;
    end        
end

P_B = max(P_pcc/etaOW, P_pcc*etaOW);

%% ECONOMIC EVALUATION

% Economic data
c_el = 1.8;                         % [NTD/kWh]
dh_capacity_price = 550;            % [NTD/MW/h]
performance_price = 350;            % [NTD/MW/h]
net_exchange = sum(P_pcc)*dt;       % [MWh]

revenue = Pserv*T*dt*(dh_capacity_price + performance_price) + min(net_exchange,0) * c_el * 1000; % [NTD]

rev_MW = revenue * 365 / (T*dt/24) / Pserv * 0.03;           % [€]
cycl_num = sum(abs(P_B*dt)) / Emax * 365 / (T*dt/24) / 2;    % [-]
fprintf('Battery cycles: %f\n', cycl_num)

%% PLOTS
% 2 main subplots: power&energy, frequency
figure
subplot(2,1,1)
hold on
plot(Tst,P_pcc)
yyaxis right 
plot(Tst,E(1:end-1)/Emax, 'r')
hold on
plot(Tst,ELref(1:end-1),'--k')
plot(Tst,ELup(1:end-1),'--g')
plot(Tst,ELdwn(1:end-1),'--y')
legend('Power [MW]','Energy [MWh]')
title('Dispatch profile')
subplot(2,1,2)
plot(Tst,f)
title('Frequency')
linkaxes(get(gcf,'children'),'x')

% Plot sorted power
figure(2)
plot(T*dt-(1:T)*dt, (sort(abs(P_B))))
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

scatter(f(1:100:end),P_pcc(1:100:end),5,E(1:100:end-1)/Emax,'filled')
colorbar

xlabel('Frequency [Hz]')
ylabel('% nominal power [-]')

% %%
% output=timetable(Tst',P_pcc_target_AFR',P_pcc',P_B');
% 
% % Tstnew=datetime(2018,12,01,0,0,0):seconds(4):datetime(2018,12,31,23,59,59);
% 
% % output=retime(output,Tstnew,'mean');
% 
% output=output(isbetween(Tst,datetime(2018,12,13,0,0,0),datetime(2018,12,20,23,59,59))',:);

