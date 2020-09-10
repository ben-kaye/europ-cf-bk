function rotors = init_graphic(offset, radii, N_sides)    
    beta = 0:2*pi/N_sides:2*pi;
    p = radii*[ cos(beta); sin(beta) ]; %
    
    rotors = [ p + offset; p + [ offset; -offset ]; p - offset; p + [ -offset; offset ] ];
end

