function plot_Cali_map(lonmin,lonmax,latmin,latmax)

load data/faults.mat;
load data/CaliCoast.dat;

% All of California
if (nargin < 4)
    lonmin = -125;
    lonmax = -114;
    latmin = 32;
    latmax = 43;
end

gcf;
coast = plot(CaliCoast(:,1),CaliCoast(:,2));
set(coast,'Color',[0 0 0],'linewidth',3)
hold on;

flts = plot(faults(:,1),faults(:,2));
set(flts,'color',[0 0 0])
axis([lonmin,lonmax,latmin,latmax])
axis square
