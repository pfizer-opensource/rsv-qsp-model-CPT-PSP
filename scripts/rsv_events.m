function [position,isterminal,direction] = rsv_events(t,x,VL_peak_PBO,event_floor)
    if (~exist('event_floor', 'var'))
        event_floor = 1;
    end

    position = VL_peak_PBO*x(4)-event_floor; % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = -1;   % locates only zeros where the event function is decreasing
end