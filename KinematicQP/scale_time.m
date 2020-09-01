function tau = scale_time(t, dist2)
    a = 0.5; % scale factor in x
    offset = 1;
    if dist2 <= offset
        tau = t;
    else
        tau = t/(1 + 1/a*(dist2 - offset)^2);
    end    
end

