%% Input Launch Parameters

% Launch site

launchsite = [53.104684, -2.924442]; % The Griffin Inn, Rossett
launchsite_err = 1*2*10^-4; % in degrees: 1*30-meter radius overall error
launch_altitude = 12.8; % launch altitude in m
launch_altitude_err = 2.3; % launch altitude error in m

% Vehicle Parameters

total_mass = 2705; % total mass (g)
total_mass_err = 23.9; % error (g)
balloon_mass_with_excess = 1292; % mass of balloon (g)
balloon_mass_with_excess_err = 3; % error(g)
payload_mass = 1505; % mass of everything & the balloon's excess mass from the theoretical value (g)
payload_mass_err = 23.9; % payload mass with excess error (g)
payload_mass_without_excess = 1413; % mass of everything minus balloon (g)
payload_mass_without_excess_err = 23.7; % error (g)

% Measurement Uncertainties

neck_err = 150; % uncertainty in the neck lift measurement at launch (g)

% Descent Rate: the balloon remains are modelled as exerting exactly a 0 
% force, with an effective mass of that of the payload

load('parachuteall_lim_3.5_11.mat') % load parachute speed distribution data
minspeed = 3.5; % in m/s
maxspeed = 11; % in m/s


% Launch Time - make sure it's not in the past or more than 1 week in thefuture!

launch_date = '27-Oct-2018'; % 'YYYY-Mon-DD' - Need to specify month with 3-letter name -
launch_time = '13:00'; % 'HH:MM'
launch_time_error = 10; % error in minutes


%% SIMULATION PARAMETERS

nsim = 50; % number of simulations (typically ~100 will give a reasonable 
% answer, 1000K+ will give very accurate stats - use ~300 for forecasting and
% ~2000 for final launch predictions)