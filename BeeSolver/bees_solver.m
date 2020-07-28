
% 6 dof in Q and 3 dof in R
% Define Qw = [ Q11 Q33 Q44 Q66 Q77 Q99 ]; Rw = [ R11 R22 R44];
Qwlower = [ 1e-1 1e-1 1e-2 1e-2 1e-2 1e-1 ];
Qwupper = [ 5e2 5e2 5e1 5e1 5e2 30 ];

Rwlower = [ 1e-1 1e-1 1e-1];
Rwupper = [ 2e2 2e2 2e2 ];

% parameters
n = 20; % no. scouts
m = 10; % no. patches
e = 4; % no. elite patches

nep = 10; % bees at elite patches
nsp = 7; % bees at other patches

ngh = 0.1; % radius of patch 
kmax = 30; % max iters

patches = {};
scouts = {};


best_cost = inf;

for g = 1:n
    scouts{g} = new_scout(Qwlower, Qwupper,Rwlower, Rwupper);    
end

function bee = new_bee(patch)
    r = patch.Radius;
    c = patch.Center;
    
    Qw = (c - r) + 2*r*rand(1, 6);
    Rw = (c - r) + 2*r*rand(1, 3);
    
    [Qa, Ra] = get_QR_actual(Qw,Rw);
    cost = get_fitness(Qa,Ra);
    
    bee = struct('Qw',Qw,'Rw',Rw, 'Cost',cost);
end

function scout = new_scout(Qwlower, Qwupper, Rwlower, Rwupper)
    Qw = Qwlower + (Qwupper-Qwlower).*rand(1,6);
    Rw = Rwlower + (Rwupper-Rwlower).*rand(1,3);
    [Qa, Ra] = get_QR_actual(Qw,Rw);
    cost = get_fitness(Qa,Ra);
    scout = struct('Qw',Qw,'Rw',Rw, 'Cost', cost);
end

function patch = make_patch(scout, nbees, radius)
    center = scout.Qw;
    radius = ngh;
    bees = {};
    best_cost = scout.cost;
    
    patch = struct('Center', center, 'Radius', radius, 'Bees', bees, 'Cost', best_cost);
end

function [Qactual, Ractual] = get_QR_actual(Qw, Rw)
    Qactual = [ Qw(1), Qw(1), Qw(2), Qw(3), Qw(3), Qw(4), Qw(5), Qw(5), Qw(6) ];
    Ractual = [ Rw(1), Rw(2), Rw(2), Rw(3) ];
end