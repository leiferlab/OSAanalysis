% reads osb data directly from serial port, play realtime movie

%% test: talk to teensy

% (run with rui_serialCommunicationTest.ino)

teensy = serial('COM3');
fopen(teensy);

% expected data structure: 2 uint32 numbers = 8 bytes total
ncol = 2;
nbytes = ncol*4;
counter = 0;

fwrite(teensy,'\n'); % send begin command

while counter<10
    if teensy.BytesAvailable>=nbytes
        a = fread(teensy,ncol,'uint32');
        if teensy.ValuesReceived==ncol
            t0 = a(1);
        end
        disp([sprintf('%0.0f',a(1)) 's, rand=' sprintf('%0.0f',a(2))]);
        counter = counter+1;
    end
end

fclose(teensy);
clear('teensy');

%% prepare .txt file path to save data to

% dstdir = 'G:\My Drive\Gershow Lab\temp_data';
% addpath(genpath(dstdir));

dstdir = pwd;

fname = [datestr(datetime,'yyyymmddHHMMSS') '_osb.txt']; % osa: array; osb: bar
[fname,dstdir] = uiputfile('.txt','Save data to',fullfile(dstdir,fname));

%% visualizer: set up axes and bar objects

fig = figure('Units','pixels','Position',[100 100 1450 800]);
% OS CO2/TVOC readings: 1st column
A(1) = axes('Units','pixels','Position',[50,600,400,150]);
B(1) = bar(A(1),[(1:8);(9:16)]'); title(A(1),'CO2 (ppm)');
% A(2) = axes('Units','pixels','Position',[50,425,400,150]);
B(2) = bar(A(1),(8:1)');
linkaxes(A(1:2),'y');
A(3) = axes('Units','pixels','Position',[50,225,400,150]);
B(3) = bar(A(3),NaN(8,1)); title(A(3),'TVOC (ppb)');
A(4) = axes('Units','pixels','Position',[50,50,400,150]);
B(4) = bar(A(4),NaN(8,1));
linkaxes(A(3:4),'y');
% OS H2/Ethanol raw signals: 2nd column
A(5) = axes('Units','pixels','Position',[500,600,400,150]);
B(5) = bar(A(5),NaN(8,1)); title(A(5),'raw H2');
A(6) = axes('Units','pixels','Position',[500,425,400,150]);
B(6) = bar(A(6),NaN(8,1));
linkaxes(A(5:6),'y');
A(7) = axes('Units','pixels','Position',[500,225,400,150]);
B(7) = bar(A(7),NaN(8,1)); title(A(7),'raw ethanol');
A(8) = axes('Units','pixels','Position',[500,50,400,150]);
B(8) = bar(A(8),NaN(8,1));
linkaxes(A(7:8),'y');
% HS T/RH readings: 3rd column
A(9) = axes('Units','pixels','Position',[1000,600,400,150]);
B(9) = bar(A(9),NaN(8,1)); title(A(9),'T (C)'); A(9).YLim = [0 100];
A(10) = axes('Units','pixels','Position',[1000,225,400,150]);
B(10) = bar(A(10),NaN(8,1)); title(A(10),'RH (%)'); A(10).YLim = [0 100];

%% TEMP load offline data

[fname,fdir] = uigetfile(fullfile(dstdir,'*.txt'),"select osa.txt");
D = importdata(fullfile(fdir,fname));
D = D.data;
t0 = D(1,1);

%% visualizer: TEMP offline update with saved data

for i=1:size(D,1)
    if(~ishandle(fig))
        break;
    end
    databuf = D(i,:);
    updatePlot(A,B,databuf,t0);
%     if i<size(D,1)
%         dt = (D(i+1,1)-D(i,1))/1000;
%         pause(dt);
%     end
%     fprintf(fid,'%u\t%u\t%u\t%u\t%u\t%u\t%u\r\n',databuf);
end

% fclose(fid);

%% visualizer: online update once every sensor, saving data

% establish serial communication
teensy = serial('COM3');
teensy.InputBufferSize = 4096; % max. inbound data size per measurement cycle = 2240B
fopen(teensy);

% expected data structure
ncol = 7;
nbytes = ncol*4; % uint32 type

fwrite(teensy,'\n'); % send begin command
fid = fopen(fullfile(dstdir,fname),'w'); % create data file
fprintf(fid,'ms_timer\tcode\tsensor_index\tsensor_ID\tdata_1\tdata_2\tvalid\r\n');

while ishandle(fig)
    if teensy.BytesAvailable>=nbytes
        databuf = fread(teensy,ncol,'uint32');
        if teensy.ValuesReceived==ncol
            t0 = databuf(1);
        end
        updatePlot(A,B,databuf,t0);
        fprintf(fid,'%u\t%u\t%u\t%u\t%u\t%u\t%u\r\n',databuf);
    end
end

fclose(fid);
fclose(teensy);
delete(teensy);
clear fig A B teensy;

%% serial object control

teensy = serial('COM3');
% teensy.InputBufferSize = 4096; % max. inbound data size per measurement cycle = 2240B
fopen(teensy);

fclose(teensy);
delete(teensy);
clear teensy;

%% STEP 1 visualizer: axis limit control

% set either/both limit to Inf/NaN to auto update ylim instead

ylim_co2 = [0 Inf]; % datasheet: 400~60000 ppm
ylim_tvoc = [0 Inf]; % datasheet: 0~60000 ppb
ylim_h2 = [10000 20000]; % datasheet: 0~65535 (uint16 range)
ylim_ethanol = [10000 20000]; % datasheet: 0~65535 (uint16 range)
ylim_t = [0 100]; % datasheet: -40~100
ylim_rh = [0 100]; % datasheet: 0~100

ylims = [ylim_co2;ylim_tvoc;ylim_h2;ylim_ethanol;ylim_t;ylim_rh];

%% STEP 2 visualizer: online update once every measurement cycle, saving data

dstdir = pwd;
fname = [datestr(datetime,'yyyymmddHHMMSS') '_osb.txt']; % osa: array; osb: bar
[fname,dstdir] = uiputfile('.txt','Save data to',fullfile(dstdir,fname));

[fig,A,B] = prepareFigureWindow('ylims',ylims);

% begin serial communication
teensy = serial('COM4');
fopen(teensy);

% expected data structure
ncol = 7;
nbytes = ncol*4; % uint32 type

% power-on sequence
disp('power-on sequence...');
fwrite(teensy,'p'); tic;
while toc<=5
    if teensy.BytesAvailable>=1
        err = fread(teensy,1,'uint8'); % 0=ok
        assert(err==0,'ERROR: failed to power-on all sensors! Please retry');
        break;
    end
end
if teensy.BytesAvailable>0
    fread(teensy,teensy.BytesAvailable); % flush input buffer
end

% OS init sequence
disp('OS init sequence...');
fwrite(teensy,'i'); tic;
while toc<=5
    if teensy.BytesAvailable>=1
        err = fread(teensy,1,'uint8'); % 0=ok
        if err~=0
            warning('unexpected output from OS init sequence!');
        end
        % ATTENTION: this is only indicating that the init sequence has
        % finished; we don't really have a way to know whether it's
        % successful or not
        break;
    end
end
if teensy.BytesAvailable>0
    fread(teensy,teensy.BytesAvailable); % flush input buffer
end

% create .txt file to save sensor data
fid = fopen(fullfile(dstdir,fname),'w');
fprintf(fid,'ms_timer\tcode\tsensor_index\tsensor_ID\tdata_1\tdata_2\tvalid\r\n');

% flush input buffer, prepare to receive sensor data
if teensy.BytesAvailable>0
    fread(teensy,teensy.BytesAvailable);
end

% measurement sequence
disp('begin measurements...');
fwrite(teensy,'m');
while ishandle(fig)
    if teensy.BytesAvailable>=nbytes
        d = fread(teensy,ncol,'uint32');
        if teensy.ValuesReceived<ncol*2 % TODO need to be smarter here
            t0 = d(1);
        end
        switch d(2)
            case 1 % OS measurements
                dd = [d fread(teensy,[ncol 15],'uint32')]'; % load the next 15 rows
                badpts = ~dd(:,end);
                ydata = dd(:,5:6);
                ydata(badpts,:) = NaN;
                ydata = reshape(ydata,[8 4])';
                ta = 1:2; tb = 1:4; % set target plots to update
            case 2 % OS raw signals
                dd = [d fread(teensy,[ncol 15],'uint32')]';
                badpts = ~dd(:,end);
                ydata = dd(:,5:6);
                ydata(badpts,:) = NaN;
                ydata = reshape(ydata,[8 4])';
                ta = 3:4; tb = 5:8;
            case 3 % HS measurements
                dd = [d fread(teensy,[ncol 7],'uint32')]';
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
        for j=length(ta)
            if any(~isfinite(ylims(ta(j),:)))
                A(ta(j)).YLimMode = 'auto';
            end
        end
        et = (dd(end,1)-t0)/1000;
        xlabel(A(5),['elapsed time: ' sprintf('%0.3f',et) 's'],...
            'FontSize',24,'Position',[4.5 -50]);
        if et>15
            xlabel(A(6),'15s initialization complete',...
                'FontSize',24,'Position',[4.5 -50]);
        end
        drawnow;
        for j=1:size(dd,1)
            fprintf(fid,'%u\t%u\t%u\t%u\t%u\t%u\t%u\r\n',dd(j,:));
        end
    end
end

% close data file
fclose(fid);
clear fig A B t0;

% random char to break measurement cycle? otherwise might have trouble
% beginning reset sequence
fwrite(teensy,'a'); 
pause(2);

% reset sequence
if teensy.BytesAvailable>0
    fread(teensy,teensy.BytesAvailable); % flush input buffer
end
disp('reset sequence...');
fwrite(teensy,'r'); tic;
while toc<=5
    if teensy.BytesAvailable>=1
        err = fread(teensy,1,'uint8'); % 0=ok
        if err~=0
            warning('unexpected output from reset sequence!');
        end
        break;
    end
end

% end serial communication
fclose(teensy);
delete(teensy);
clear teensy;

%% test: send char command to teensy

teensy = serial('COM3');
fopen(teensy);

fwrite(teensy,'a');
pause(0.1);
while teensy.BytesAvailable>0
    b = fread(teensy,1,'uint8');
    disp(['teensy state: ' num2str(b)]);
end

fclose(teensy);
delete(teensy);
clear teensy;

%% test: control flow

tic;
while toc<=5
    toc;
    if toc>=2
        tmp = toc;
    end
    assert(tmp<=3);
    return;
end