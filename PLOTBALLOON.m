
%% Plot Preliminary Results
close all;
figure(1);


lat_all = cell([1, length(ascent_lat)]);
lon_all = cell([1, length(ascent_lat)]);
alt_all = cell([1, length(ascent_lat)]);

for i = 1:length(ascent_lat)
    
    lat_all{i} = [ascent_lat{i}, descent_lat{i}];
    lon_all{i} = [ascent_lon{i}, descent_lon{i}];
    alt_all{i} = [ascent_alt{i}, descent_alt{i}];
    
end

% DEFINE SCALE
xminlon = min(cell2mat(lon_all))-0.05;
xmaxlon = max(cell2mat(lon_all))+0.05;
ymaxlat = max(cell2mat(lat_all))+0.05;
yminlat = min(cell2mat(lat_all))-0.05;
ymeanlat = mean(cell2mat(lat_all));


shortestindex = min([100, length(ascent_lon)]);

for plotn = 1:shortestindex
    
    plot(ascent_lon{plotn},ascent_lat{plotn},'r-')
    hold on;
    plot(descent_lon{plotn},descent_lat{plotn},'b-')
    plot(landing_lon(plotn),landing_lat(plotn),'kx')
    legend('Ascent', 'Descent', 'Landing','Location','NorthWest');
    titlestr = sprintf('Flight Prediction [%s]', launch_time_str);
    title(titlestr);
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    
end



%% Fit Error Ellipse to Landing Site

h1 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.682689492137086);
set(h1,'LineWidth',1);
set(h1,'Color','g');
set(h1,'LineStyle','--');
h2 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9545);
set(h2,'LineWidth',1);
set(h2,'Color','g');
set(h2,'LineStyle','--');
h3 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9973);
set(h3,'LineWidth',1);
set(h3,'Color','g');
set(h3,'LineStyle','--');
xlim([xminlon, xmaxlon]);
ylim([yminlat, ymaxlat]);
plot_google_map('Scale',2,'Alpha',0.6,'MapType','hybrid', 'Refresh', 0,'AutoAxis',0)
daspect([(1/cos(deg2rad(ymeanlat))) 1 1])
hold off;

figure(2)


if ymaxlat-yminlat < xmaxlon-xminlon
    
    meanxlon = (xminlon+xmaxlon)/2;
    rangexlon = xmaxlon-xminlon;
    ymaxlat = ymeanlat+(9/16)*(1/2)*(xmaxlon-xminlon);
    yminlat = ymeanlat-(9/16)*(1/2)*(xmaxlon-xminlon);
    xminlon = meanxlon-rangexlon/2*(1/cos(deg2rad(ymeanlat)));
    xmaxlon = meanxlon+rangexlon/2*(1/cos(deg2rad(ymeanlat)));
    
else
    
    meanlat = (yminlat+ymaxlat)/2;
    rangeylat = ymaxlat-yminlat;
    xmeanlon = mean(cell2mat(lon_all));
    xmaxlon = xmeanlon+(16/(9*1))*(1/2)*rangeylat*(1/cos(deg2rad(ymeanlat)));
    xminlon = xmeanlon-(16/(9*1))*(1/2)*rangeylat*(1/cos(deg2rad(ymeanlat)));

end

plot(landing_lon,landing_lat,'kx')
axis equal tight;
lonplotrange = [mean(landing_lon)-3*std(landing_lon) mean(landing_lon)+3*std(landing_lon)];
latplotrange = [mean(landing_lat)-3*std(landing_lat) mean(landing_lat)+3*std(landing_lat)];
xlim([xminlon, xmaxlon]);
ylim([yminlat, ymaxlat]);

plot_google_map('Scale',2,'Resize',2,'Alpha',0.5,'MapType','hybrid','Refresh', 0,'AutoAxis',0);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
daspect([(1/cos(deg2rad(ymeanlat))) 1 1])


%% AMAZING STUFF!!

figure(3)

h1 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.682689492137086);
handle12 = gcf;

axesObjs = get(handle12, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata1 = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata1 = get(dataObjs, 'YData');

close()

figure(3)

h2 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9545);
handle12 = gcf;

axesObjs = get(handle12, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata2 = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata2 = get(dataObjs, 'YData');

close()

figure(3)
h3 = error_ellipse(cov(landing_lon,landing_lat),[mean(landing_lon),mean(landing_lat)],'conf',0.9973);
handle12 = gcf;

axesObjs = get(handle12, 'Children');  %axes handles     
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
xdata3 = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata3 = get(dataObjs, 'YData');
close()

figure(3)


% DEFINE SCALE
xminlon = min(landing_lon)-0.01;
xmaxlon = max(landing_lon)+0.01;
ymaxlat = max(landing_lat)+0.01;
yminlat = min(landing_lat)-0.01;
ymeanlat = mean(landing_lat);

if ymaxlat-yminlat < xmaxlon-xminlon
    
    meanxlon = (xminlon+xmaxlon)/2;
    rangexlon = xmaxlon-xminlon;
    ymaxlat = ymeanlat+(9/16)*(1/2)*(xmaxlon-xminlon);
    yminlat = ymeanlat-(9/16)*(1/2)*(xmaxlon-xminlon);
    xminlon = meanxlon-rangexlon/2*(1/cos(deg2rad(ymeanlat)));
    xmaxlon = meanxlon+rangexlon/2*(1/cos(deg2rad(ymeanlat)));
    
else
    
    meanlat = (yminlat+ymaxlat)/2;
    rangeylat = ymaxlat-yminlat;
    xmeanlon = mean(cell2mat(lon_all));
    xmaxlon = xmeanlon+(16/(9*1))*(1/2)*rangeylat*(1/cos(deg2rad(ymeanlat)));
    xminlon = xmeanlon-(16/(9*1))*(1/2)*rangeylat*(1/cos(deg2rad(ymeanlat)));

end

% plot intensity distribution
nsimres = length(landing_lon);

if nsimres < 400
    scalefactor = 1;
elseif and(nsimres >= 400, nsimres < 1000)
    scalefactor = 2;
elseif nsimres >= 1000;
    scalefactor = 6;
elseif nsimres >= 6000;
    scalefactor = 9;
end

[xlength,ylength] = rat(range(landing_lon)/range(landing_lat),3e-4);
orgrid = (xlength*ylength)*6*scalefactor; % make the scalefactor 1 for small numbers of samples < 500
xres = round(xlength*nsimres/orgrid)*15;
yres = round(ylength*nsimres/orgrid)*15;

if xres < 3
    xres = 10;
end

if yres < 3
    yres = 10;
end
    
resolution = [xres yres]; % set number of horizontal and vertical bins (pixels)
intensitydata = [transpose(landing_lon) transpose(landing_lat)];
temperaturemap = hist3(intensitydata,resolution);
Vq = interp2(temperaturemap, 3)'; % 3D histogram the x and y interaction data
Vq = Vq/max(Vq(:)); % normalise
xb = linspace(min(intensitydata(:,1)),max(intensitydata(:,1)),size(temperaturemap,1)+1); % create the x and y axes
yb = linspace(min(intensitydata(:,2)),max(intensitydata(:,2)),size(temperaturemap,1)+1); 
% build a colourmap image of the intensity distribution
% the log function extends the range of interactions far
% beyond the centre of the distribution
alphamask = Vq>0.1;
hnd = imagesc(xb,yb,Vq);
axis equal;
xlim([xminlon, xmaxlon]);
ylim([yminlat, ymaxlat]);
set(hnd, 'AlphaData', alphamask);
%set(hnd, 'FaceAlpha', 0.1);
%alpha(hnd,.3);
titlestr = sprintf('Landing Probability Distribution [%s]', launch_time_str);
title(titlestr);
ylabel('Latitude (deg)');
hold on;
xlabel('Longitude (deg)');
abitofahack = colorbar();
ylabel(abitofahack, 'Relative Landing Probability','FontSize',12);
colormap('default');
plot_google_map('Scale',2,'Alpha',1,'MapType','hybrid','AutoAxis',0);


hh1 = plot(xdata1, ydata1, 'g--', 'LineWidth', 1);
hh2 = plot(xdata2, ydata2, 'y--', 'LineWidth', 1);
hh3 = plot(xdata3, ydata3, 'r--', 'LineWidth', 1);

legend([hh1 hh2 hh3],{'1\sigma','2\sigma','3\sigma'})
daspect([(1/cos(deg2rad(ymeanlat))) 1 1])

