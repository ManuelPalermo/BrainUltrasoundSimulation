clear; clc;

% create the computational grid
Nx = 180;           % number of grid points in the x (row) direction
Ny = 180;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;    % [m/s]

% define the array of time points [s]
t_end = 20e-6;      % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);



%%%%%%%%%%%%%%%% sensors
n_points = 30;

%step=5; sensor_points = [1:step:step*n_points]+10000;
sensor_points = sort(randi(prod([Nx,Ny])/2,[n_points,1]));

vector_mask = zeros(prod([Nx,Ny]),1);
vector_mask(sensor_points)=1;

sensor.mask = reshape(vector_mask, [Nx,Ny]);




%%%%%%%%%%%%%%%%%%%%% source
source.p_mask = zeros(Nx, Ny);
source.p_mask(round(Ny/1.2), round(Nx/1.5)) = 1;

sampling_freq = 1/kgrid.dt;     % [Hz]
tone_burst_freq = 1e6;          % [Hz]
tone_burst_cycles = 2;

source.p = 10*toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);


%% Calculate sensor delays to concentrate wave on given position

input_args = {'PlotLayout', true, ... 
              'PlotPML', false, ...
              'DisplayMask', source.p_mask + sensor.mask,...
              'RecordMovie', true, 'MovieName', 'PingTarget',...
              'DataCast', 'single'};
          
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


% plot the time series recorded at each array element
figure;
stackedPlot(kgrid.t_array * 1e6, sensor_data);
xlabel('Time [\mus]');
ylabel('Time Series Recorded At Each Element');

[max_val, max_pos] = max(sensor_data,[],2);
%travel_time = max_pos*kgrid.dt;
travel_time = max_pos;
sensor_delays = max(travel_time)-travel_time;


%% Simulate concentrated wave on a given position

% swap sensor to source and monitor pressure on all grid
focus_point_mask = source.p_mask;
source.p_mask = sensor.mask;
sensor.mask = ones(Nx,Ny);
sensor.record = {'p', 'p_max'};

% add delay to source points
source.p = 1*toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
                        'SignalOffset', sensor_delays);

% create new simulation
input_args = {'PlotLayout', true, ... 
              'PlotPML', false, ...
              'DisplayMask', source.p_mask + focus_point_mask,...
              'RecordMovie', true, 'MovieName', 'FocusTarget',...
              'DataCast', 'single'};
          
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


figure;
p_max = reshape(sensor_data.p_max, Nx, Ny);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p_max);
h = colorbar;
xlabel(h, '[MPa]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Max Acoustic Pressure Amplitude');




