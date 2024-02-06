function[ascent_rate, burst_altitude, neck_lift_sim] = ASCENTRATE(launch_vol, balloon_mass_with_excess_sim, payload_mass_without_excess_sim, neck_err)
%% ASCENTRATE
% Ascent rate, bust altitude and neck lift calculator.
%
% INPUT: 
%   launch_vol: pre-set launch volume to reach 5 m/s at mean performance (m^3)
%   balloon_mass_with_excess_sim: simulated total balloon mass (g)
%   payload_mass_without_excess_sim: simulated total payload + cord + parachute mass
%   neck_err: error in the measurement of neck lift at launch, limited by
%   instruments to ~5 g typically.
%
% OUTPUT: 
%   launch_vol: pre-set launch volume to reach 5 m/s at mean performance 
%   balloon_mass_with_excess_sim: simulated total balloon mass
%   payload_mass_without_excess_sim: simulated total payload + cord + parachute mass
%   neck_err: error in the measurement of neck lift at launch, limited by
%   instruments to ~5 g typically.

% Gas & Atmospheric Parameters

gas_dens = 0.1786; % Helium density(Kg/m^3)
air_dens = 1.2050; % Air density at 0C, 101 kPa (Kg/m^3)
atm_model = 7238.3; % Atmospheric model scale height (m)
grav_acc = 9.80665; % Mean gravitational acceleration (m/s^2)

% Balloon Launch Parameters

launch_dia = 2*power((3*launch_vol)/(4*pi),(1/3)); % Launch diameter (m)
launch_area = pi*power((launch_dia/2),2); % Launch area

% Simulated Burst Diameter

a = 0.3112; % Burst diameter uses dia/m = a*(x/grams)^b	where a_err = 0.05165; and b_err = 0.0226;
b = 0.4658;

burst_dia = a*balloon_mass_with_excess_sim^b; % in m
burst_dia_err = 0.1238*2; % 2 sigma error from model RMSE (m) THIS MAY BE A SLIGHT UNDERESTIMATION
burst_dia_sim = normrnd(burst_dia,burst_dia_err); % simulated burst diameter

% Simulated Burst Volume

burst_vol = (4/3)*pi*power(burst_dia_sim/2,3); % simulated burst volume
burst_volume_ratio = burst_vol/launch_vol; % simulated burst to launch volume ratio

% Simulated Burst Altitude

burst_altitude = -(atm_model*log(1/burst_volume_ratio)); % simulated burst altitude in meters

% Simulated Lift

gross_lift = launch_vol*(air_dens-gas_dens)*1000; % in grams
free_lift = gross_lift - balloon_mass_with_excess_sim - payload_mass_without_excess_sim; % in grams

neck_lift = free_lift+payload_mass_without_excess_sim; % in grams
neck_lift_sim = normrnd(neck_lift,neck_err); % simulated lift due to precision of luggage scale in grams

free_lift_sim = neck_lift_sim-payload_mass_without_excess_sim; % simulated free lift

free_lift_force = grav_acc*free_lift_sim/1000; % simulated free lift force in N 

% Simulated Drag Coefficient - insuficient data here!

cd = 0.25;
cd_err = 0.005; % again, the error is unknown, but mostly in the range 0.25-0.30
cdsim = normrnd(cd,cd_err);

% Simulated Ascent Rate

ascent_rate = sqrt(free_lift_force/(0.5*cdsim*air_dens*launch_area)); % m/s

% Other interesting equations:
% time_to_burst = (burst_altitude/ascent_rate)/60; % min

end