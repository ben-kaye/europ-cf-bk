function rotors = init_graphic()
    offset = 0.2;
    radii = 0.1;
    
    
    N_sides = 13;
    beta = 0:2*pi/N_sides:2*pi;
    p = radii*[ cos(beta); sin(beta) ];
    
    rotors = [ p + offset; p + [ offset; -offset ]; p - offset; p + [ -offset; offset ] ];
end

