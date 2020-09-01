function tau = scale_time(t, dist2, scale, offset)
%     if dist2 <= offset
%         tau = t;
%     else
%         tau = t/(1 + 1/scale*(dist2 - offset)^2);
%     end    
    tau = t/(1 + scale/dist2);
end

