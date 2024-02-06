function[vol_guess] = LAUNCHVOLUME(desired_speed, mean_balloon_mass_with_excess, mean_payload_mass_without_excess)
%% LAUNCHVOLUME
% Launch volume calculator, using a hill-climbing loop.
%
% INPUT: 
%   desired_speed: usually 5 m/s
%   mean_balloon_mass_with_excess: mean total balloon mass (g)
%   mean_payload_mass_without_excess: mean total payload + cord + parachute mass (g)
%
% OUTPUT: 
%   vol_guess: pre-set launch volume to reach 5 m/s at mean performance (m^3)

% Gas & Atmospheric Parameters

gas_dens = 0.1786; % Helium density(Kg/m^3)
air_dens = 1.2050; % Air density at 0C, 101 kPa (Kg/m^3)
grav_acc = 9.80665; % Mean gravitational acceleration (m/s^2)
cd = 0.25; % Mean Drag Coefficient
vol_guess = 4; % Initial Volume Guess
step = 0.001; % Precision



speed_diff = 1;

while abs(speed_diff) > step
    
    % correct value 
    if speed_diff > 0 % if speed is higher 
        vol_guess = vol_guess - 0.001; % inflate balloon less
    elseif speed_diff < 0
        vol_guess = vol_guess + 0.001; % inflate balloon more
    end
    
    launch_dia = 2*power((3*vol_guess)/(4*pi),(1/3)); % Launch diameter (m)
    launch_area = pi*power((launch_dia/2),2); % Launch area
    gross_lift = vol_guess*(air_dens-gas_dens)*1000; % in grams
    free_lift = gross_lift - mean_balloon_mass_with_excess - mean_payload_mass_without_excess; % in grams
    neck_lift = free_lift+mean_payload_mass_without_excess; % in grams
    free_lift_force = grav_acc*free_lift/1000; % mean free lift force in N 
    
    
    calculated_speed = sqrt(free_lift_force/(0.5*cd*air_dens*launch_area)); % m/s
    
    speed_diff = calculated_speed - desired_speed;
    
    
end


% fprintf('\nPreset Speed = %0.5f m/s\nNeck Lift = %0.1f g\n', calculated_speed, neck_lift)