function ctrls_1 = sat_ctrls(ctrls, min_ctrl, max_ctrl)
    ctrls1 = max(ctrls, min_ctrl);
    ctrls1 = min(ctrls1, max_ctrl);
end

