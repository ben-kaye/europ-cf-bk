function r = ref_lookup(r, t, step_sz, path_index)
    switch path_index
        case 1
            if t < 0.6
                r(4) = r(4) - 1.1*step_sz;
            else
                if t < 1
                    r(4) = -0.9;
                    r(6) = -0.05;
                else 
                    r(4) = 0;
                    r(6) = 0;
                end
            end
        case 0
            r(4) = 0.5;
            r(6) = 0;
    end
end

