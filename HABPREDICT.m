%% HABPREDICT V1.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Roke Cepeda-Arroita    %
%  University of Manchester  %
%        January 2017        %
%       V1.2 October 2018    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This program queries the CUSF Landing Predictor 2.5 to calculate landing
% error ellipses of 1-sigma (68% confidence), 2-sigma (95%) and 3-sigma
% (99.7%) regions. This is done by simulating the following initial
% parameters: balloon burst altitude, neck lift uncertainties, vehicle mass
% uncertainties, the terminal parachute descent speed distribution and the
% time and site of the launch. It relies on four external scripts: RANDPDF,
% LAUNCHVOLUME, ASCENTRATE and error_ellipse, as well as a file containing 
% the probability distribution of the terminal descent speeds. The program 
% will output a launch card with useful information such as the amount of 
% helium needed, the distance to the landing site and the magnitude of the 
% search area, as well as csv files containing the landing ellipses that 
% can be made into a map using Google "My Maps" for the recovery process 
% and planning of the chase route. The API can return about 7 simulations 
% every second, but 200 should be enough to get a good idea of the landing 
% area. DISCLAIMER: this doesn't take into account parachute glide (which 
% can lead to errors of up to ~1 km), or early balloon bursts. GFS weather
% estimates can also be biased above ~30 km. Simulation results can be 
% visualised using the script PLOTBALLOON.

% ***V1.2 Log***
% Added config file with all the variables, so change that and when you are
% happy with your variables run HABPREDICT. Then PLOTBALLOON for visuals.


%% Clear Workspace

clearvars;
close all;


%% Import Configuration file

CONFIG % This runs CONFIG.m and copies all configuration options over here


%% Simulate Parameters

h = waitbar(0,'Simulating parameters...');

% create descent speeds

descent_rates = RANDPDF(meandist,index,[1,nsim*20]);
descent_rates = descent_rates(and(descent_rates>3.5,descent_rates<11)); % apply cutoff
mean_descent_rate_all = mean(descent_rates);
std_descent_rate_all = std(descent_rates);
descent_rates = descent_rates(1:nsim);


for i = 1:nsim
    
    %disp(i)
    
    waitbar(i/nsim);

    % Simulated Balloon & Payload Masses

    balloon_mass_with_excess_sim = normrnd(balloon_mass_with_excess,balloon_mass_with_excess_err);
    payload_mass_without_excess_sim = normrnd(payload_mass_without_excess,payload_mass_without_excess_err);

    % Simulated Launch Time

    launch_time_str = datetime(sprintf('%s %s:00', launch_date, launch_time));
    launch_time_str2 = datestr(launch_time_str,'dd.mmm.yyyy.HH.MM');
    launch_time_delay = round(normrnd(0,3*60)); % round to nearest second
    launch_time_sim = launch_time_str + seconds(launch_time_delay);
    launch_time_unchanged = datevec(launch_time_str);
    launch_time_sim = datevec(launch_time_sim);

    % Simulated Descent Rate

    descent_rate_sim = descent_rates(i);

    % Simulated Launch Altitude

    launch_altitude_sim = normrnd(launch_altitude,launch_altitude_err); % launch altitude in m

    % Simulated Launch Site

    launchsite_sim(1) = normrnd(launchsite(1),launchsite_err);
    launchsite_sim(2) = normrnd(launchsite(2),launchsite_err);


    %% Calculate Ascent Speed

    % Calculate Mean Launch Volume to Achieve 4.8 m/s

    [launch_vol] = LAUNCHVOLUME(4.8, balloon_mass_with_excess, payload_mass_without_excess);

    % Flight Parameters by Fixing the Launch Volume and Simulating Everything Else

    [ascent_rate_sim(i), burst_altitude_sim(i), neck_lift_sim(i)] = ASCENTRATE(launch_vol, balloon_mass_with_excess_sim, payload_mass_without_excess_sim, neck_err);

    %% Create Query String

    % Correct Launch Site Longitude (0,360)

    if launchsite_sim(2) < 0
        launchsite_sim(2) = launchsite_sim(2)+360;
    end

    % Create Variable Strings

    latstr = sprintf('launch_latitude=%0.5f', launchsite_sim(1));
    lonstr = sprintf('launch_longitude=%0.5f', launchsite_sim(2));
    altstr = sprintf('launch_altitude=%0.1f', launch_altitude_sim);
    datetimestr = sprintf('launch_datetime=%.4d-%.2d-%.2dT%.2d%%3A%.2d%%3A%.2d%%2B00:00', launch_time_sim);
    ascentstr = sprintf('ascent_rate=%0.3f', ascent_rate_sim(i));
    burstaltstr = sprintf('burst_altitude=%0.0f', burst_altitude_sim(i));
    descentstr = sprintf('descent_rate=%0.3f', descent_rate_sim);

    % Create URL
    

   	APIurl = sprintf('http://predict.cusf.co.uk/api/v1/?%s&%s&%s&%s&%s&%s&%s',latstr,lonstr,altstr,datetimestr,ascentstr,burstaltstr,descentstr);


    
    % EXAMPLE: http://predict.cusf.co.uk/api/v1/?launch_latitude=50.0&launch_longitude=0&launch_altitude=72&launch_datetime=2016-11-29T18%3A32%3A01%2B00:00&ascent_rate=5.1&burst_altitude=33000&descent_rate=5.3


    %% Request API and Load Data
    
    keep_trying = 1;

    while keep_trying == 1
        keep_trying = 0;
             try
                simdat = webread(APIurl); % copy API data
             catch
             	disp('An error occurred while retrieving information from the internet.');
             	disp('Execution will continue.');
             	keep_trying = 1;
             end
             
    end
     
    % Correct Longitudes to Avoid Jumps

    ascent_lon{i} = [simdat.prediction(1).trajectory.longitude];
    ascent_lat{i} = [simdat.prediction(1).trajectory.latitude];
    ascent_alt{i} = [simdat.prediction(1).trajectory.altitude];
    descent_lon{i} = [simdat.prediction(2).trajectory.longitude];
    descent_lat{i} = [simdat.prediction(2).trajectory.latitude];
    descent_alt{i} = [simdat.prediction(2).trajectory.altitude];
    alltimes = [simdat.prediction(1).trajectory.datetime simdat.prediction(2).trajectory.datetime];
    flight_times = strsplit(alltimes,'Z');
    flight_times = flight_times(1:length(flight_times)-1);
    start_time = str2double(strsplit(flight_times{1,1},{'-','T',':'}));
    for counterz = 1:length(flight_times)
        all_time_stamps = str2double(strsplit(flight_times{counterz},{'-','T',':'}))-start_time; 
        timevectorp(counterz) = (all_time_stamps(6)+60*all_time_stamps(5)+60*60*all_time_stamps(4)+60*60*24*all_time_stamps(3))/60; % in minutes
    end
    all_time{i} = timevectorp;
    end_time = str2double(strsplit(flight_times{1,length(flight_times)-1},{'-','T',':'}));
    pred_time1 = [simdat.request.dataset];
    pred_time(i) = datenum(str2double(strsplit(pred_time1(1:length(pred_time1)-1),{'-','T',':'})));
    total_flight_time(i) = etime(end_time,start_time);
    
    ascent_lon{i}(ascent_lon{i} > 180) = ascent_lon{i}(ascent_lon{i} > 180)-360;
    descent_lon{i}(descent_lon{i} > 180) = descent_lon{i}(descent_lon{i} > 180)-360;
    
    landing_lon(i) = descent_lon{i}(length(descent_lon{i}));
    landing_lat(i) = descent_lat{i}(length(descent_lat{i}));
    

end

close(h);



%% THROW OUT BAD VALUES, REWRITE SIMULATION DATA

max_displacement = 0.1;
badpredictions = [];

% Identify rogue paths

for i = 1:nsim

    if any(abs(diff(ascent_lon{i}))>max_displacement) == 1
        
        badpredictions = [badpredictions, i];
        
    elseif any(abs(diff(descent_lon{i}))>max_displacement) == 1
        
        badpredictions = [badpredictions, i];
        
    end
    
    
end

goorpredictions = setdiff(1:nsim, badpredictions);

% Rewrite data

landing_lon = landing_lon(goorpredictions);
landing_lat = landing_lat(goorpredictions);
total_flight_time = total_flight_time(goorpredictions);

ascent_lon(badpredictions) = [];
ascent_lat(badpredictions) = [];
descent_lon(badpredictions) = [];
descent_lat(badpredictions) = [];


%% GET MEAN PATH

% Create Query String

% Correct Launch Site Longitude (0,360)

[ascent_rate_mean, burst_altitude_mean, neck_lift_sim_mean] = ASCENTRATE(launch_vol, balloon_mass_with_excess_sim, payload_mass_without_excess_sim, 0);

launchsitemod = launchsite;
if launchsite(2) < 0
    launchsitemod(2) = launchsite(2)+360;
end

% Create Variable Strings

latstr = sprintf('launch_latitude=%0.5f', launchsitemod(1));
lonstr = sprintf('launch_longitude=%0.5f', launchsitemod(2));
altstr = sprintf('launch_altitude=%0.1f', launch_altitude);
datetimestr = sprintf('launch_datetime=%.4d-%.2d-%.2dT%.2d%%3A%.2d%%3A%.2d%%2B00:00', launch_time_unchanged);
ascentstr = sprintf('ascent_rate=%0.3f', mean(ascent_rate_sim));
burstaltstr = sprintf('burst_altitude=%0.0f', mean(burst_altitude_sim));
descentstr = sprintf('descent_rate=%0.3f',mean_descent_rate_all);

% Create URL

APIurl = sprintf('http://predict.cusf.co.uk/api/v1/?%s&%s&%s&%s&%s&%s&%s',latstr,lonstr,altstr,datetimestr,ascentstr,burstaltstr,descentstr);

% EXAMPLE: http://predict.cusf.co.uk/api/v1/?launch_latitude=50.0&launch_longitude=0&launch_altitude=72&launch_datetime=2016-11-29T18%3A32%3A01%2B00:00&ascent_rate=5.1&burst_altitude=33000&descent_rate=5.3


% Request API and Load Data

simdat = webread(APIurl); % copy API data

% Correct Longitudes to Avoid Jumps


ascent_altitude = [simdat.prediction(1).trajectory.altitude];
ascent_lon_mean = [simdat.prediction(1).trajectory.longitude];
ascent_lat_mean = [simdat.prediction(1).trajectory.latitude];
descent_altitude = [simdat.prediction(2).trajectory.altitude];
descent_lon_mean = [simdat.prediction(2).trajectory.longitude];
descent_lat_mean = [simdat.prediction(2).trajectory.latitude];

ascent_lon_mean(ascent_lon_mean > 180) = ascent_lon_mean(ascent_lon_mean > 180)-360;
descent_lon_mean(descent_lon_mean > 180) = descent_lon_mean(descent_lon_mean > 180)-360;

meanpath = zeros(length(ascent_lon_mean)+length(descent_lon_mean),2);
meanpath(:,1) = [ascent_lat_mean,descent_lat_mean];
meanpath(:,2) = [ascent_lon_mean,descent_lon_mean];
meanaltitude = [ascent_altitude,descent_altitude];

% Work out total travelled distance

totaldistance = 0;
horizontalspeed = zeros([1,length(meanaltitude)-1]);
travelled_km = zeros([1,length(meanaltitude)-1]);

for i = 1:length(meanaltitude)-1

    horizontaldistance = 111.120*distance(meanpath(i,1),meanpath(i,2),meanpath(i+1,1),meanpath(i+1,2)); % in km
    horizontalspeed(i) = horizontaldistance/60*1000;
    horzvertdistance = sqrt(horizontaldistance^2+(meanaltitude(i+1)/1000-meanaltitude(i)/1000)^2);
    
    totaldistance = totaldistance + horzvertdistance;
    travelled_km(i) = totaldistance;
    
end

horizontalspeed = [horizontalspeed,0];

% Work out total subtended angle in the ascent path of the balloon

totangle = 0;

for i = 1:length(meanpath)-length(descent_lat_mean)-2
    
    a = [meanpath(i+1,:) 0]-[meanpath(i,:) 0]; % initial vector
    b = [meanpath(i+2,:) 0]-[meanpath(i+1,:) 0]; % subsequent vector
    
    newangle(i) = rad2deg(atan2(norm(cross(a,b)),dot(a,b)));
    
    totangle = totangle + newangle(i);
    
end

possible_angle_factor = std(diff(newangle));

totangle2 = 0;

for i = length(ascent_lat_mean)+1:length(meanpath)-2
    
    a = [meanpath(i+1,:) 0]-[meanpath(i,:) 0]; % initial vector
    b = [meanpath(i+2,:) 0]-[meanpath(i+1,:) 0]; % subsequent vector
    
    newangle = rad2deg(atan2(norm(cross(a,b)),dot(a,b)));
    
    totangle2 = totangle2 + newangle;
    
end

%% SIMULATE EARLY BURSTS

early_delta_h = linspace(0,-20000,8);
early_burst_altitudes = mean(burst_altitude_sim)+early_delta_h(2:length(early_delta_h));

for early_i = 1:length(early_burst_altitudes)

% Create Variable Strings

latstr = sprintf('launch_latitude=%0.5f', launchsitemod(1));
lonstr = sprintf('launch_longitude=%0.5f', launchsitemod(2));
altstr = sprintf('launch_altitude=%0.1f', launch_altitude);
datetimestr = sprintf('launch_datetime=%.4d-%.2d-%.2dT%.2d%%3A%.2d%%3A%.2d%%2B00:00', launch_time_unchanged);
ascentstr = sprintf('ascent_rate=%0.3f', mean(ascent_rate_sim));
burstaltstr = sprintf('burst_altitude=%0.0f', early_burst_altitudes(early_i));
descentstr = sprintf('descent_rate=%0.3f',mean_descent_rate_all);

% Create URL

APIurl = sprintf('http://predict.cusf.co.uk/api/v1/?%s&%s&%s&%s&%s&%s&%s',latstr,lonstr,altstr,datetimestr,ascentstr,burstaltstr,descentstr);

% EXAMPLE: http://predict.cusf.co.uk/api/v1/?launch_latitude=50.0&launch_longitude=0&launch_altitude=72&launch_datetime=2016-11-29T18%3A32%3A01%2B00:00&ascent_rate=5.1&burst_altitude=33000&descent_rate=5.3


% Request API and Load Data

simdat = webread(APIurl); % copy API data

% Correct Longitudes to Avoid Jumps


descent_lon_early = [simdat.prediction(2).trajectory.longitude];
descent_lat_early = [simdat.prediction(2).trajectory.latitude];

descent_lon_early(descent_lon_early > 180) = descent_lon_early(descent_lon_early > 180)-360;

early_lon(early_i) = descent_lon_early(length(descent_lon_early));
early_lat(early_i) = descent_lat_early(length(descent_lat_early));



end;



earlybursts = zeros(length(early_lon),3);
earlybursts(:,1) = early_lat;
earlybursts(:,2) = early_lon;
earlybursts(:,3) = early_delta_h(2:length(early_delta_h));


%% GET HANDLE ON ERROR ELLIPSE
% SEMIMAJOR, MINOR AXES, ANGLE, CENTRE!

landing_centre = [mean(landing_lat),mean(landing_lon);mean(landing_lat),mean(landing_lon)];

figure(1)
h1 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.682689492137086);
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
close()

sigma1_x = xdata(1:3:length(xdata));
sigma1_y = ydata(1:3:length(ydata));

figure(2)
h2 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9545);
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
close()

sigma2_x = xdata(1:3:length(xdata));
sigma2_y = ydata(1:3:length(ydata));

figure(3)
h3 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9973);
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
close()

sigma3_x = xdata(1:3:length(xdata));
sigma3_y = ydata(1:3:length(ydata));

%% Save as CSV

csvstr = sprintf('%sMEANPATH.csv',launch_time_str2);
csvwrite(csvstr,meanpath);

csvstr = sprintf('%sEARLYBURST.csv',launch_time_str2);
csvwrite(csvstr,earlybursts);

csvstr = sprintf('%sCENTRE.csv',launch_time_str2);
csvwrite(csvstr,landing_centre);

sigma_1 = zeros(length(sigma1_x),2);
sigma_1(:,2) = sigma1_x;
sigma_1(:,1) = sigma1_y;

csvstr = sprintf('%sSIGMA1.csv',launch_time_str2);
csvwrite(csvstr,sigma_1);

sigma_2 = zeros(length(sigma2_x),2);
sigma_2(:,2) = sigma2_x;
sigma_2(:,1) = sigma2_y;

csvstr = sprintf('%sSIGMA2.csv',launch_time_str2);
csvwrite(csvstr,sigma_2);

sigma_3 = zeros(length(sigma3_x),2);
sigma_3(:,2) = sigma3_x;
sigma_3(:,1) = sigma3_y;

csvstr = sprintf('%sSIGMA3.csv',launch_time_str2);
csvwrite(csvstr,sigma_3);

areaoflanding = areaint(sigma_1(:,1),sigma_1(:,2))*5.100644719097883e+08; % in km^2



%% Print Test Results

fprintf('\nINPUT ***********************')
fprintf('\n%s ± %0.0f min\n', launch_time_str,launch_time_error);
fprintf('(%0.6f,%0.6f) ± %0.0f m\n', launchsite, 111120*distance(launchsite(1),launchsite(2),launchsite(1)+launchsite_err/sqrt(2),launchsite(2)+launchsite_err/sqrt(2)) );
fprintf('LAUNCH ALTITUDE = %0.0f ± %0.0f m\n', launch_altitude, launch_altitude_err);
fprintf('BALLOON MASS = %0.1f ± %0.1f g\n', balloon_mass_with_excess,balloon_mass_with_excess_err);
fprintf('PAYLOAD MASS = %0.1f ± %0.1f g\n', payload_mass_without_excess,payload_mass_without_excess_err);
fprintf('DESCENT RATE = %0.2f ± %0.2f m/s\n', mean_descent_rate_all,std_descent_rate_all);
fprintf('\nOUTPUT ***********************')
fprintf('\nDATASET = %s (±%0.0f) \n', datestr(datevec(mean(pred_time)), 'dd-mmm HH:MM PM') , std(pred_time)*24*60); 
fprintf('MEAN HELIUM VOLUME = %0.3f m^3\n', launch_vol);
fprintf('NECK LIFT = %0.2f ± %0.2f g\n', mean(neck_lift_sim),std(neck_lift_sim));
fprintf('ASCENT RATE = %0.2f ± %0.2f m/s\n', mean(ascent_rate_sim),std(ascent_rate_sim));
fprintf('BURST ALTITUDE = %0.0f ± %0.0f m\n', mean(burst_altitude_sim),std(burst_altitude_sim));
fprintf('FLIGHT DURATION = %s ± %0.0f min\n', datestr(mean(total_flight_time)/(3600*24), 'HH:MM') ,std(total_flight_time)/60);
fprintf('AVE LANDING ACCURACY = %0.1f km\n', 111.120*distance(mean(landing_lat),mean(landing_lon),mean(landing_lat)+std(landing_lat),mean(landing_lon)+std(landing_lon)) );
fprintf('LAT LANDING ACCURACY = %0.1f km\n', 111.120*distance(mean(landing_lat),mean(landing_lon),mean(landing_lat)+std(landing_lat),mean(landing_lon) ));
fprintf('LON LANDING ACCURACY = %0.1f km\n', 111.120*distance(mean(landing_lat),mean(landing_lon),mean(landing_lat),mean(landing_lon)+std(landing_lon)) );
fprintf('1-SIGMA SEARCH AREA = %0.1f km^2\n', areaoflanding);
fprintf('DISTANCE FROM LAUNCH = %0.0f km\n', 111.120*distance(mean(landing_lat),mean(landing_lon),launchsite(1),launchsite(2)));
fprintf('TRAVELLED BY BALLOON = %0.0f km\n', totaldistance);
fprintf('ASCENT HEADING SHIFT = %0.0f deg\n', totangle);
fprintf('DESCENT HEADING SHIFT = %0.0f deg\n\n', totangle2);