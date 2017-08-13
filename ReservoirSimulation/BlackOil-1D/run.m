addpath(genpath(pwd));

%% This is one time step
%  With more steps, update old values

%% Environment
% phase ID
OIL = 1;
WAT = 2;

% direction 
L = 1;
R = 2;

% grid
Ncell = 3;
Nface = Ncell-1;

% grid layout
% center: 1 - 2 - 3 - 4 - 5
% face:   - 1 - 2 - 3 - 4 -

DX = 20;    % [ m ]
DY = 20;    % [ m ]
DZ = 20;    % [ m ]
DT = 100; % [ sec ]

DV = DX * DY * DZ; % [ m3 ]
%% Constant Properties
constant_water_density = 1000; %[ kg/m3 ]
Co = 0.03; %[ kg/(m3.Pa) ]
porosity = 0.3;
perm = 5E-12; %[ m2 ]
viscosity( OIL ) = 1.1E-03; %[ Pa.s ]
viscosity( WAT ) = 0.8E-03; %[ Pa.s ]
flowcoeff = FlowCoeff( OIL, WAT, perm, viscosity, DX );

%% Wellbore
WI = Well_Index( DX, DZ, perm );

% Produce with constant BHP
% wid_prd  = [ 3 ];
% BHP_PRD = [ 20000 ]; % [ Pa ]

wid_prd  = [  ];
BHP_PRD = [  ]; % [ Pa ]

% Inject with constant rate
wid_inj = [ 1 ];
const_Qw_inj = [ -80 ];

% wid_inj = [  ];
% const_Qw_inj = [  ];

%% Storage
incr = length( wid_inj );
P  = zeros( Ncell + incr, 1 );
Sw = zeros( Ncell, 1 );

% P_old  = zeros( Ncell, 1 );
% Sw_old = zeros( Ncell, 1 );

Res = zeros( 2*Ncell+incr, 1 );
J = zeros( 2*Ncell+incr, 2*Ncell+incr );

%% Initialization
Pi = 2.5E+04; %[ Pa ]
Swi = 0.5;
for i = 1 : Ncell
    P(i)  = Pi + (i-1) * 1000;
    Sw(i) = Swi + (i-1) * 0.02;
end

for i = 1 : length( wid_inj )
    P( Ncell + incr, 1 ) = 25000;
end

% Center & Interfacial Properties (including derivatives)
[dens, rperm, dens_DP, rperm_DS] = Update_Center( Ncell, OIL, WAT, P, Sw, Co, constant_water_density );
[DP,densf, rpermf] = Update_Interface( Nface, OIL, WAT, P, dens, rperm );

P_old = P;
Sw_old = Sw;

print_flag = 0; % no print
%% Residual Calculation
Res = Compute_Res( Res, Ncell, Nface, OIL, WAT, P, Sw, Sw_old, P_old, dens, densf, rperm, rpermf, DP, flowcoeff, Co, constant_water_density, porosity, viscosity, DT, WI, wid_inj, const_Qw_inj, wid_prd, BHP_PRD, DV, print_flag );

%% Jacobian
J = Compute_Jacobian( J, Ncell, Nface, OIL, WAT, L, R, P, Sw, dens, densf, dens_DP, rperm, rpermf, rperm_DS, DP, flowcoeff, porosity, viscosity, DT, WI, wid_inj, wid_prd, BHP_PRD, DV );

%% Newton
ITER = 4;
iter = 0;

TOL = 1E-03;
norm_Res = 1;
norm_Del = 1;

while( (norm_Res >= TOL || norm_Del >= TOL) && iter < ITER )   
    iter = iter + 1;
    
    % solve
    delta = -J \ Res;
    
    norm_Del = norm( delta );
    NORM_DEL( iter, 1 ) = norm_Del;
    
    norm_Res = norm( Res );
    NORM_RES( iter, 1 ) = norm_Res;
    
    if( iter == 4 )
        pr(Res);
    end
    
    % Appleyard chop
    for j =1 : Ncell
        delSw = delta( 2*j, 1 );
        if abs( delSw ) > 0.2
            delta( 2*j, 1 ) = sign( delSw ) * 0.2;
        end
    end
    
% This is the bug~~!!!!!!
% Old values DO NOT CHANGE within newton step~~!!!!!!
%     P_old = P;
%     Sw_old = Sw;
    
    % update Pressure
    for i = 1 : length(P)
        P(i) = P(i) + delta( 2*i-1, 1 );
    end

    % update Saturation
    for i = 1 : length(Sw)
        Sw(i) = Sw(i) + delta( 2*i, 1 );
        if( Sw(i) <= 0 )
            Sw(i) = 0.001;
        elseif Sw(i) >= 1
            Sw(i) = 0.999;
        end
    end   
    
    [dens, rperm, dens_DP, rperm_DS] = Update_Center( Ncell, OIL, WAT, P, Sw, Co, constant_water_density );
    [DP,densf, rpermf] = Update_Interface( Nface, OIL, WAT, P, dens, rperm );
    
    if iter == -1
        print_flag = 1;
    else
        print_flag = 0;
    end
    
    Res = Compute_Res( Res, Ncell, Nface, OIL, WAT, P, Sw, Sw_old, P_old, dens, densf, rperm, rpermf, DP, flowcoeff, Co, constant_water_density, porosity, viscosity, DT, WI, wid_inj, const_Qw_inj, wid_prd, BHP_PRD, DV, print_flag );
    J = Compute_Jacobian( J, Ncell, Nface, OIL, WAT, L, R, P, Sw, dens, densf, dens_DP, rperm, rpermf, rperm_DS, DP, flowcoeff, porosity, viscosity, DT, WI, wid_inj, wid_prd, BHP_PRD, DV );
 
end



