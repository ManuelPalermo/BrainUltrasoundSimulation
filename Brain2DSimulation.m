

%% Carregar modelo cerebro e objetivo a focar
clear; clc;

% carregar modelo cerebro
%brain_model = niftiread('brain_model.nii');
brain_model = load('brain_model_skull.mat'); brain_model = brain_model.brain_model;
brain_model = brain_model(25:225, 35:274, 39:239); % dimunuir espaço inutil

% definir uma slice
model = squeeze(brain_model(85,:,:));

% visualizar o modelo
figure;
imshow(model);

%% Definir o meio e suas propriedades

% create the 2D - computational grid
[Nx, Ny] = size(model);   % number of grid points in the X/Y direction
dx = 1e-3;                % grid point spacing in the X direction [m]
dy = 1e-3;                % grid point spacing in the Y direction [m]
dz = 1e-3;                % grid point spacing in the Z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);


% define the properties of the propagation medium
%tissue      speed[m/s]  density(Kg/m3)  absortion[dB/(MHz.cm)] img-pixel-values
%air          - 343.0      -  1.20        - 0.0004              - [0]*
%water        - 1475.0     -  1000        - 0.05?               - [0]*
%midbrain     - 1546.3     -  1000        - 0.6                 - [21-39] U [51-78]
%white matter - 1552.5     -  1050        - 0.6                 - [40-50]
%grey matter  - 1500.0     -  1100        - 0.6                 - [81-220]
%cSpinalFluid - 1475.0     -  1000        - 0.05                - [1-9]
%scalp        - 1540.0     -  1000        - 0.1                 - [10-20]
%skull        - 3476.0     -  1979        - 2.7                 - [221-255]
%* caso seja usado agua ou ar como meio exterior

% speed [m/s]
medium.sound_speed = 1500 * ones(Nx, Ny);               % default
%medium.sound_speed(model==0)=343.0 ;                    % air
medium.sound_speed(model==0)=1504.0;                    % water
medium.sound_speed(model>=21 & model<=78)=1546.3;       % midbrain
medium.sound_speed(model>=40 & model<=50)=1552.5;       % white matter
medium.sound_speed(model>=81 & model<=220)=1500.0;      % grey matter
medium.sound_speed(model>=1  & model<=9)=1475.0;        % cerebroSpinalFluid
medium.sound_speed(model>=10 & model<=20)=1540.0;       % scalp
medium.sound_speed(model>=221)=3476.0;                  % skull
% density [Kg/m3]
medium.density = 1000 * ones(Nx, Ny);                   % default
%medium.density(model==0)=1.20;                          % air
medium.density(model==0)=1000;                          % water
medium.density(model>=21 & model<=78)=1075;             % midbrain
medium.density(model>=40 & model<=50)=1050;             % white matter
medium.density(model>=81 & model<=220)=1100;            % grey matter
medium.density(model>=1  & model<=9)=1000.0;            % cerebroSpinalFluid
medium.density(model>=10 & model<=20)=1000.0;           % scalp
medium.density(model>=221)=1969.0;                      % skull
% absortion [dB/(MHz^y cm)]
medium.alpha_power = 1.5;                               % default
medium.alpha_coeff = 0.75 * ones(Nx, Ny);               % default
%medium.alpha_coeff(model==0)=1.6;                       % air
medium.alpha_coeff(model==0)=0.05;                      % water
medium.alpha_coeff(model>=21 & model<=78)=0.6;          % midbrain
medium.alpha_coeff(model>=40 & model<=50)=0.6;          % white matter
medium.alpha_coeff(model>=81 & model<=220)=0.6;         % grey matter
medium.alpha_coeff(model>=1  & model<=9)=0.05;          % cerebroSpinalFluid
medium.alpha_coeff(model>=10 & model<=20)=0.1;          % scalp
medium.alpha_coeff(model>=221)=2.7;                     % skull

% tempo e step de simulação
kgrid.makeTime(medium.sound_speed);


%% Define ultrassound source and sensors

n_elements = 20;                        % grid points - impar
source.p_mask = zeros(Nx, Ny);

source2use = 3;
%Source Position:
    % 1 - Array de pontos sequenciais colados ao topo do cerebro
    % 2 - Pontos aleatorios no lado de cima do cerebro
    % 3 - Pontos completamente aleatorios

if source2use==1
    %%%%%%% Probe com varios elementos colados em cima do cerebro %%%%%%%%%
    array_step=6;
    array_center_pos_y = size(model,1)/2;
    for ypos=0:array_step:array_step*(n_elements-1)/2    % puts all aray elements in contact with the brain
        % lado positivo em relacao ao centro
        xpos = find(model(array_center_pos_y+ypos,:)>0,1,'last');
        source.p_mask(array_center_pos_y+ypos, xpos) = 1;
        % lado negativo em relacao ao centro
        xpos = find(model(array_center_pos_y-ypos,:)>0,1,'last');
        source.p_mask(array_center_pos_y-ypos, xpos) = 1;
    end
end

if source2use==2
    %%%%%%%%%%%% Fontes aleatorias do lado de cima do cerebro %%%%%%%%%%%%%
    source_points = sort(randi(numel(model)/3,[n_elements,1])) + numel(model)*(2/3);
    vector_mask = zeros(prod([Nx,Ny]),1);
    vector_mask(source_points)=1;
    source.p_mask = reshape(vector_mask, [Nx,Ny]);
end

if source2use==3
    %%%%%%%%%%%% Fontes aleatorias na matriz %%%%%%%%%%%%%
    source_points = sort(randi(numel(model),[n_elements,1]));
    vector_mask = zeros(prod([Nx,Ny]),1);
    vector_mask(source_points)=1;
    source.p_mask = reshape(vector_mask, [Nx,Ny]);
end

ping_pressure = 20;                     % [Pa]
signal_freq = 0.25e6;                   % [Hz]
ping_burst_cycles = 1;
source.p = ping_pressure * toneBurst(1/kgrid.dt, signal_freq, ping_burst_cycles);


% sensors
sensor.mask = ones(Nx,Ny);
sensor.record = {'p', 'p_max'};

%% Simulation

                 
input_args = {'PlotLayout', false, ... 
              'PlotPML', false, ...
              'DisplayMask', source.p_mask | model==255,...
              'RecordMovie', true, 'MovieName', 'SimulationVideo',...
              'DataCast', 'single'};
          
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

figure;
p_max = reshape(sensor_data.p_max, Nx, Ny);
%p_max(model<21 | model>221) = 0;           % ignore points outside the brain
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p_max);
h = colorbar; xlabel(h, '[Pa]');
title('Max Acoustic Pressure');

