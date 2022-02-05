clear
clc

% load('\\10.36.8.12\Rd_Sys\01_COMMESSE Q\Q1045 - PPC EMS\02_TechSaleConsulting\Hoping Projects\Data\TCC_frequency.mat')

% Profiles including the SoC reference according to "standard" SoC
% management strategy (independent composition of peak shaving & dReg)
load('ProfiliHoping')
% Profiles of frequency and peak shaving setpoint
load('ProfiliHoping_NoSoCRef')

%% OPTIONS AND PARAMETERS

OpMode = "EdReg";      %EdReg,dReg05,dReg025

Tst=datetime(2018,12,01,0,0,0):seconds(1):datetime(2018,12,31,23,59,59);
T=length(Tst);

% Storage and Service parameters
Emax = 286;         % [MWh]
Pserv = 100;        % [MW]

ELtarget=0.1;
deadband=0.05;
dt=1/3600;
etaOW=0.94;

f(:)=60;


%% 

% SoCref = nan(1,T);
% if ~strcmp(OpMode,EdReg)
%     % Killing peak shaving (only dReg)
%     SoCref(:)=0.5;
%     Pps_AC(:)=0;
% end

if strcmp(OpMode,"dReg0.25")
    % dReg 0.25 parameters
    % Map of prescribed response
    upcurve=[100 100 52 9 9 -52 -100 -100];
    lowcurve=[100 100 52 -9 -9 -52 -100 -100];
    freq=[59 59.75 59.86 59.98 60.02 60.14 60.25 61];
else
    % dReg 0.5 parameters (valid also for EdReg)
    % Map of prescribed response
    upcurve=[100 100 48 9 9 -48 -100 -100];
    lowcurve=[100 100 48 -9 -9 -48 -100 -100];
    freq=[59 59.5 59.75 59.98 60.02 60.25 60.5 61];
end

% Discharge-oriented response
dchresponse=interp1(freq,upcurve,f)/100*Pserv;
% Charge-oriented response
chresponse=interp1(freq,lowcurve,f)/100*Pserv;

% figure
% plot(dchresponse,'--k')
% hold on
% plot(chresponse,'--k')
%%
T=length(dchresponse);
E = nan(1,T+1);
ELref = nan(1,T+1);
P_soc_mgnt = zeros(1,T);
P_pcc_target_AFR = nan(1,T);
E(1) = Emax*0.5;
ELref(1) = ELtarget;
ELup = ELtarget+deadband;
ELdwn = ELtarget-deadband;

norespflag=sign(dchresponse.*chresponse)==-1;
percnoresp=sum(norespflag)/T;
norespset=inf(T,1);
norespset(norespflag)=0;
% Lowest allowed response
zeroresponse=min(abs([dchresponse';chresponse';norespset'])).*sign(60-f');

% Pps_DC=max(Pps_AC/etaOW,Pps_AC*etaOW);


hflag = 0;
P_pcc = nan(1,T);

% figure
% histogram(zeroresponse,-100.5:1:100.5,'normalization','probability')
% sum(abs(zeroresponse/3600),'omitnan')/250

%% Dispatch Simulation

hflag=0;

for t=1:T
    
    if strcmp(OpMode,"EdReg")
        % Adjusting reference level for AFR SoC management based on the
        % integration of the peak shaving setpoint
        ELtarget = ELref(t);
        ELup = ELtarget+deadband;
        ELdwn = ELtarget-deadband;
    end

    flaghist(t)=hflag;
    
    Pmin = -(Emax - E(t))/etaOW/dt;
    Pmax = E(t)*etaOW/dt;
    
    if hflag == 0
        
        if Pps_AC(t) == 0
            P_pcc_target_AFR(t) = zeroresponse(t);
        elseif Pps_AC(t) > 0
            P_pcc_target_AFR(t) = chresponse(t);  
            P_soc_mgnt(t) = chresponse(t)-zeroresponse(t);
        elseif Pps_AC(t) < 0
            P_pcc_target_AFR(t) = dchresponse(t);
            P_soc_mgnt(t) = dchresponse(t)-zeroresponse(t);
        end            
        
    elseif hflag == 1
        
        P_pcc_target_AFR(t) = dchresponse(t);
        
    elseif hflag == -1
        
        P_pcc_target_AFR(t) = chresponse(t);
        
    end
    
    Pps_AC_adj(t) = Pps_AC(t) + P_soc_mgnt(t);
    Pps_DC_adj(t) = max(Pps_AC_adj(t)/etaOW,Pps_AC_adj(t)*etaOW);

    Pps_DC(t) = max(Pps_AC(t)/etaOW,Pps_AC(t)*etaOW);
    
    P_pcc_target(t) = P_pcc_target_AFR(t) + Pps_AC(t);
    
    P_pcc(t) = min([max([P_pcc_target(t),Pmin,-Pserv]),Pmax,Pserv]);
    
    SBSPM(t) = 100 - abs(P_pcc(t)-P_pcc_target(t))/Pserv*100;
    
    E(t+1) = E(t) - max(P_pcc(t)/etaOW,P_pcc(t)*etaOW)*dt;
    ELref(t+1) = ELref(t) - Pps_DC_adj(t)*dt/Emax;
%     ELref(t+1) = ELref(t) - Pps_DC(t)*dt/Emax;

    if (E(t+1)/Emax)>ELup
        hflag=1;
    elseif E(t+1)/Emax<ELdwn
        hflag=-1;
    elseif ((E(t+1)/Emax)<ELtarget&&hflag==1)||((E(t+1)/Emax)>ELtarget&&hflag==-1)
        hflag=0;
    end        
    
    
end

P_B = max(P_pcc/etaOW,P_pcc*etaOW);

HVAC_cons=sum((max(P_B*(1-0.96),0)-min(P_B*(1-0.96)/0.96,0))*dt)/3.5;

%%
figure
subplot(2,1,1)
plot(Tst,P_pcc)
hold on
% plot(E)
yyaxis right 
plot(Tst,E(1:end-1)/Emax)
hold on
plot([Tst(1),Tst(end)],[ELup ELup],'--k')
plot([Tst(1),Tst(end)],[ELdwn ELdwn],'--k')
legend('Power [MW]','Energy [MWh]')
title('Dispatch profile')
% plot(flaghist)
subplot(2,1,2)
plot(Tst,f)
title('Frequency')
linkaxes(get(gcf,'children'),'x')

%% Economic evaluation

c_el = 1.8;                         % [NTD/kWh]
dh_capacity_price = 550;            % [NTD/MW/h]
performance_price = 350;            % [NTD/MW/h]
net_exchange = sum(P_pcc)*dt;       % [MWh]

revenue = Pserv*T*dt*(dh_capacity_price + performance_price) + min(net_exchange,0) * c_el * 1000; % [NTD]

rev_MW = revenue*365/(T*dt/24)/Pserv*0.03           % [€]
cycl_num = sum(abs(P_B*dt))/Emax*365/(T*dt/24)/2    % [-]

figure
plot((sort(abs(P_B))),T*dt-[1:T]*dt)
xlabel('Battery power [MW]')
ylabel('# of operating hours at higher battery power [h]')

figure
plot(freq(2:end-1),upcurve(2:end-1),'--k','linewidth',2)
hold on
plot(freq(2:end-1),lowcurve(2:end-1),'--k','linewidth',2)

plot(freq(2:end-1),lowcurve(2:end-1)+50,'--k','linewidth',2)
plot(freq(2:end-1),upcurve(2:end-1)+50,'--k','linewidth',2)

plot(freq(2:end-1),lowcurve(2:end-1)-50,'--k','linewidth',2)
plot(freq(2:end-1),upcurve(2:end-1)-50,'--k','linewidth',2)

scatter(f(1:100:end),P_pcc(1:100:end),5,E(1:100:end-1)/Emax,'filled')
colorbar


%%
output=timetable(Tst',P_pcc_target_AFR',P_pcc',P_B');

% Tstnew=datetime(2018,12,01,0,0,0):seconds(4):datetime(2018,12,31,23,59,59);

% output=retime(output,Tstnew,'mean');

output=output(isbetween(Tst,datetime(2018,12,13,0,0,0),datetime(2018,12,20,23,59,59))',:);

