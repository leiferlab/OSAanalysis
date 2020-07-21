% read osb data from disk, play movie

% trying to use bar.YDataSource + refreshdata: nope, even slower

%% STEP 1 import data to matlab

% dstdir = 'G:\My Drive\Gershow Lab\temp_data';
dstdir = pwd;
[fname,fdir] = uigetfile(fullfile(dstdir,'*.txt'),'select osa.txt');
D = importdata(fullfile(fdir,fname));
D = D.data;

%% STEP 2 visualizer: control variables

ylim_co2 = [0 10000]; % datasheet: 400~60000 ppm
ylim_tvoc = [0 1000]; % datasheet: 0~60000 ppb
ylim_h2 = [10000 15000]; % datasheet: 0~65535 (uint16 range)
ylim_ethanol = [15000 20000]; % datasheet: 0~65535 (uint16 range)
ylim_t = [0 100]; % datasheet: -40~100
ylim_rh = [0 100]; % datasheet: 0~100

ylims = [ylim_co2;ylim_tvoc;ylim_h2;ylim_ethanol;ylim_t;ylim_rh];

%% visualizer: set up axes and bar objects and data sources

fig = figure('Units','pixels','Position',[100 100 1000 900]);

% OS readings
% CO2
A(1) = axes('Units','pixels','Position',[50,700,400,150]);
B(1:2) = bar(A(1),NaN(8,2));
A(1).Title.String = 'CO2 (ppm)';
% TVOC
A(2) = axes('Units','pixels','Position',[50,500,400,150]);
B(3:4) = bar(A(2),NaN(8,2));
A(2).Title.String = 'TVOC (ppb)';
% raw H2
A(3) = axes('Units','pixels','Position',[50,250,400,150]);
B(5:6) = bar(A(3),NaN(8,2));
A(3).Title.String = 'raw H2'; A(3).YScale = 'log';
% raw Ethanol
A(4) = axes('Units','pixels','Position',[50,50,400,150]);
B(7:8) = bar(A(4),NaN(8,2));
A(4).Title.String = 'raw ethanol'; A(4).YScale = 'log';

% HS readings
% T
A(5) = axes('Units','pixels','Position',[550,700,400,150]);
B(9) = bar(A(5),NaN(8,1));
A(5).Title.String = 'T (C)';
% RH
A(6) = axes('Units','pixels','Position',[550,250,400,150]);
B(10) = bar(A(6),NaN(8,1));
A(6).Title.String = 'RH (%)';

for i=1:length(A)
    if all(isfinite(ylims(i,:)))
        A(i).YLim = ylims(i,:);
    end
end

%% try YDataSource + refreshdata: works!

DS = NaN(10,8);
for i=1:10
    B(i).YDataSource = ['DS(' num2str(i) ',:)'];
end

targetplot = 3;
targetind = 7;
DS(targetplot,targetind) = 500;
refreshdata(A(targetplot));
ylim(A(targetplot),'auto');

%% visualizer: offline update with saved data, one sensor at a time

% NOPE - this is even slower than updatePlot.m where we just rewrite
% bar.YData each time

for i=1:size(D,1)
    % quit if figure window is closed
    if(~ishandle(fig))
        return;
    end
    % read one data entry
    databuf = D(i,:);
    % identify where to update
    targetind = databuf(3);
    % identify pair of axes to update
    switch databuf(2)
        case 1 % OS CO2/TVOC
            if targetind<=8
                targetplot1 = 1;
                targetplot2 = 3;
            else
                targetplot1 = 2;
                targetplot2 = 4;
            end
        case 2 % OS H2/Ethanol
            if targetind<=8
                targetplot1 = 5;
                targetplot2 = 7;
            else
                targetplot1 = 6;
                targetplot2 = 8;
            end
        case 3 % HS T/RH
            targetplot1 = 9;
            targetplot2 = 10;
    end
    % shift targetind to match plot
    if targetind>8
        targetind = targetind-8;
    end
    % if data entry is valid, update with value recorded, otherwise
    % update with NaN
    if databuf(end)
        val1 = databuf(5);
        val2 = databuf(6);
    else
        val1 = NaN;
        val2 = NaN;
    end
    % update linked datasource
    DS(targetplot1,targetind) = val1;
    DS(targetplot2,targetind) = val2;
    % update elapsed time and init status
    xlabel(A(9),['elapsed time: ' sprintf('%0.3f',(databuf(1)-t0)/1000) 's'],...
        'FontSize',24,'Position',[4.5 -50]);
    if ((databuf(1)-t0)/1000)>15
        xlabel(A(10),'15s initialization complete','FontSize',24,'Position',[4.5 -50]);
    end
    refreshdata(fig,'base');
    drawnow;
%     if i<size(D,1)
%         dt = (D(i+1,1)-D(i,1))/1000;
%         pause(dt);
%     end
end

%% STEP 3 visualizer: offline update with saved data, whole-board

[fig,A,B] = prepareFigureWindow('ylims',ylims);

t0 = D(1,1);
et = t0;
dloc = 0;

for i=1:size(D,1)
    if ~ishandle(fig)
        return;
    end
    if i<=dloc
        continue;
    end
    dloc = i;
    d = D(dloc,:);
    switch d(2)
        case 1 % OS measurements
            dd = D(dloc:dloc+15,:); % load the next 15 rows
            dloc = dloc+15; % update location to skip future loops
            badpts = ~dd(:,end); % mark invalid sensors
            ydata = dd(:,5:6); % extract sensor data
            ydata(badpts,:) = NaN; % take out invalid data
            ydata = reshape(ydata,[8 4])'; % reshape for convinience
            ta = 1:2; tb = 1:4; % set target plots to update
        case 2 % OS raw signals
            dd = D(dloc:dloc+15,:);
            dloc = dloc+15;
            badpts = ~dd(:,end);
            ydata = dd(:,5:6);
            ydata(badpts,:) = NaN;
            ydata = reshape(ydata,[8 4])';
            ta = 3:4; tb = 5:8;
        case 3 % HS measurements
            dd = D(dloc:dloc+7,:);
            dloc = dloc+7;
            badpts = ~dd(:,end);
            ydata = dd(:,5:6);
            ydata(badpts,:) = NaN;
            ydata = reshape(ydata,[8 2])';
            ta = 5:6; tb = 9:10;
    end
    % update bar plot data
    for j=1:length(tb)
        B(tb(j)).YData = ydata(j,:);
    end
    % update y axis limit if necessary
%     for j=length(ta)
%         if any(~isfinite(ylims(ta(j),:)))
%             A(ta(j)).YLimMode = 'auto';
%         end
%     end
    dt = (dd(end,1)-t0)/1000-et;
    et = (dd(end,1)-t0)/1000;
    xlabel(A(5),['elapsed time: ' sprintf('%0.3f',et) 's'],...
        'FontSize',24,'Position',[4.5 -50]);
    if et>15
        xlabel(A(6),'15s initialization complete',...
            'FontSize',24,'Position',[4.5 -50]);
    end
    drawnow;
    %pause(dt);
end

%% raw signal vs time

inds = find(D(:,2)==2);
h2 = D(inds,5);
ethanol = D(inds,6);
valid = D(inds,end);

% reset invalid data to NaN for plotting
h2(valid==0) = NaN;
ethanol(valid==0) = NaN;

h2 = reshape(h2,16,[]);
ethanol = reshape(ethanol,16,[]);
valid = reshape(valid,16,[]);
pt = 1:size(h2,2);

figure;
subplot(1,3,1);
plot(pt,h2(1:8,:),'b'); hold on;
plot(pt,h2(9:16,:),'r'); hold off;
ylim([10000 15000]); title('raw h2');
subplot(1,3,2);
plot(pt,ethanol(1:8,:),'b'); hold on;
plot(pt,ethanol(9:16,:),'r'); hold off;
ylim([15000 20000]); title('raw ethanol');
subplot(1,3,3);
plot(pt,valid,'.');
title('data valid');
suptitle(fname);
