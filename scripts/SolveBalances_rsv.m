function [T,X] = SolveBalances_rsv(TSTART,TSTOP,Ts,params,VL_peak_PBO, IC_PBO,drug_Finhib, drug_Ninhib,tol, event_floor)

if (~exist('event_floor', 'var'))
    event_floor = 1;
end
    
%Define RHS of model
odefun = @(t,x)rsv_model(t,x,params,drug_Finhib, drug_Ninhib);

if tol == 0
    opts = odeset('Events', @(t,x)rsv_events(t,x,VL_peak_PBO, event_floor));
else
    opts = odeset('RelTol',tol,'AbsTol',tol,'Events', @(t,x)rsv_events(t,x,VL_peak_PBO, event_floor));
end

[t_PBO_new,y_PBO_new,te,ye,ie] = ode15s(odefun,TSTART:Ts:TSTOP,...
    IC_PBO,opts);

T = [t_PBO_new];
X = [y_PBO_new];

% turn off viral compartment if VL drops below limit of detection (event-based)
if (t_PBO_new(end)<TSTOP) && ~isempty(ye)
    t1 = [];
    y1 = [];
    %Shut down viral compartment
    params(1) = 0;
    params(4) = 0;
    params(5) = 0;
    opts = odeset('Events', []);
    odefun = @(t,x)rsv_model(t,x,params,drug_Finhib, drug_Ninhib);
    
    if t_PBO_new(end) - t_PBO_new(end-1)< Ts
        TSIM = t_PBO_new(end-1):Ts:TSTOP;
        [t1,y1] = ode15s(odefun,TSIM,...
        y_PBO_new(end,:),opts);
        T = [t_PBO_new(1:end-1);t1(2:end)];
        X = [y_PBO_new(1:end-1,:);y1(2:end,:)];
    else
        TSIM = t_PBO_new(end):Ts:TSTOP;
        [t1,y1] = ode15s(odefun,TSIM,...
        y_PBO_new(end,:),opts);
        T = [t_PBO_new(1:end-1);t1(1:end)];
        X = [y_PBO_new(1:end-1,:);y1(1:end,:)];
    end
   
end

end