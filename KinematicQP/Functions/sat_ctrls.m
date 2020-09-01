function ctrls_1 = sat_ctrls(ctrls, min_ctrl, max_ctrl)
    ctrls_1 = max(ctrls, min_ctrl);
    ctrls_1 = min(ctrls_1, max_ctrl);
end

