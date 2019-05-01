

%% Carregar modelo cerebro e objetivo a focar
clear; clc;

% carregar modelo cerebro
%brain_model = niftiread('brain_model.nii');
brain_model = load('brain_model_cranio.mat'); brain_model = brain_model.brain_model;
brain_model = brain_model(25:225, 35:274, 39:239); % dimunuir espaço inutil

% definir uma slice
model = squeeze(brain_model(85,:,:));

% visualizar o modelo
figure;
imshow(model);

% definir sequencia de pontos a focar com ultrassom

%{
target_points = [[70, 80];...
                 [90, 110];...
                 [130, 90];...
                 [160, 140]]
%}
target_points = [[80, 75]];
             
% mascara com os pontos que irão ser focados
focus_points_mask = zeros(size(model));
for ipoint=1:size(target_points,1)
    focus_points_mask(target_points(ipoint,1), target_points(ipoint,2)) = 1;
end

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


%% Definir source de ondas Ultrassonoras

n_elements = 37;                        % grid points - impar
sensor.mask = zeros(Nx, Ny);

source2use = 1;
%Source Position:
    % 1 - Array de pontos sequenciais colados ao topo do cerebro
    % 2 - Array de pontos sequenciais colados ao topo e abaixo do cerebro
    % 3 - Array linear em cima do cerebro
    % 4 - Pontos aleatorios no lado de cima do cerebro
    % 5 - Pontos completamente aleatorios


if source2use==1
    %%%%%%% Probe com varios elementos colados em cima do cerebro %%%%%%%%%
    array_step=6;
    array_center_pos_y = size(model,1)/2;
    for ypos=0:array_step:array_step*(n_elements-1)/2    % puts all aray elements in contact with the brain
        % lado positivo em relacao ao centro
        xpos = find(model(array_center_pos_y+ypos,:)>0,1,'last');
        sensor.mask(array_center_pos_y+ypos, xpos) = 1;
        % lado negativo em relacao ao centro
        xpos = find(model(array_center_pos_y-ypos,:)>0,1,'last');
        sensor.mask(array_center_pos_y-ypos, xpos) = 1;
    end
end

if source2use==2
    %%%%%%% Probe com varios elementos colados em cima e baixo do cerebro %%%%%%%%%
    array_step=5;
    n_elements=(n_elements-1)/2;
    array_center_pos_y = size(model,1)/2;
    for ypos=0:array_step:array_step*(n_elements-1)/2    % puts all aray elements in contact with the brain
        % lado positivo em relacao ao centro
        xpos = find(model(array_center_pos_y+ypos,:)>0,1,'last');   % cima
        sensor.mask(array_center_pos_y+ypos, xpos) = 1;
        xpos = find(model(array_center_pos_y+ypos,:)>0,1,'first');  % baixo
        sensor.mask(array_center_pos_y+ypos, xpos) = 1;
        
        % lado negativo em relacao ao centro
        xpos = find(model(array_center_pos_y-ypos,:)>0,1,'last');   % cima
        sensor.mask(array_center_pos_y-ypos, xpos) = 1;
        xpos = find(model(array_center_pos_y-ypos,:)>0,1,'first');   % baixo
        sensor.mask(array_center_pos_y-ypos, xpos) = 1;
    end
end

if source2use==3
    %%%%%%%%%%%%%%%%%%%% Probe linear em cima do cerebro %%%%%%%%%%%%%%%%%%%%%
    array_step=1;
    array_center_pos_y = size(model,1)/2;
    max_brain_x = 0;
    for ypos=0:array_step:array_step*(n_elements-1)/2  % find outer brain position for array
        max_brain_x = max(max_brain_x, find(model(array_center_pos_y+ypos,:)>0,1,'last'));
        max_brain_x = max(max_brain_x, find(model(array_center_pos_y-ypos,:)>0,1,'last'));
    end
    sensor.mask(array_center_pos_y-ypos:array_step:array_center_pos_y+ypos, max_brain_x)=1;
end

if source2use==4
    %%%%%%%%%%%% Fontes aleatorias do lado de cima do cerebro %%%%%%%%%%%%%
    sensor_points = sort(randi(numel(model)/3,[n_elements,1])) + numel(model)*(2/3);
    vector_mask = zeros(prod([Nx,Ny]),1);
    vector_mask(sensor_points)=1;
    sensor.mask = reshape(vector_mask, [Nx,Ny]);
end

if source2use==5
    %%%%%%%%%%%% Fontes aleatorias na matriz %%%%%%%%%%%%%
    sensor_points = sort(randi(numel(model),[n_elements,1]));
    vector_mask = zeros(prod([Nx,Ny]),1);
    vector_mask(sensor_points)=1;
    sensor.mask = reshape(vector_mask, [Nx,Ny]);
end


%% Mostra todo o layout(Densidade/VelocidadeSom/Fontes/Targets)
figure;
subplot(2,2,1); % Densidade do meio
imagesc(medium.density);
colormap(gray);
b=colorbar;
title('Densidade');
xlabel(b, '[Kg/L]');

subplot(2,2,3); % Velocidade do som no meio
imagesc(medium.sound_speed);
colormap(gray);
b=colorbar;
title('Velocidade');
xlabel(b, '[m/s]');


interfaces_mask = edge(medium.density.*medium.sound_speed);
subplot(2,2,[2 4]); % Interfaces US
imshow(interfaces_mask);
colormap(gray);
title('Interfaces sigificativas US');


figure;
subplot(1,2,1); % Targets US
imshow(focus_points_mask);
colormap(gray);
title('Targets US');

subplot(1,2,2); % Fontes US
imshow(sensor.mask);
colormap(gray);
title('Fontes US');

pause(0.5);

%% Definir ping dos targets de ondas Ultrassonoras

% target ping signal
ping_pressure = 20;                     % [Pa]
signal_freq = 0.25e6;                   % [Hz]
ping_burst_cycles = 1;
source.p = ping_pressure * toneBurst(1/kgrid.dt, signal_freq, ping_burst_cycles);

input_args = {'DataCast', 'single', ...
              'PlotSim', false};
   
% calcular delays para cada ponto
for ipoint=1:size(target_points,1)
    % definir target para fazer ping
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(target_points(ipoint,1),target_points(ipoint,2))=1;
    
    % simular
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    [~, max_pos] = max(sensor_data,[],2);
    sensor_delays(ipoint,:) = max(max_pos)-max_pos;
end


%% Simulate wave focusing on multiple targets sequentialy

% swap sensor to source and monitor pressure on all grid
source.p_mask = sensor.mask;
sensor.mask = ones(Nx,Ny);
sensor.record = {'p', 'p_max'};

% definir propriedades do sinal a emitir
focus_pressure = 10;                 % [Pa]
focus_burst_cycles = 2;

total_sensor_max_pressure = zeros(Nx*Ny,1);

for ipoint=1:size(target_points,1)
    % cria sinal a emitir para focar(com delays)
    source.p = focus_pressure*toneBurst(1/kgrid.dt, signal_freq, focus_burst_cycles, ...
                                   'SignalOffset', sensor_delays(ipoint,:));
    % simular                    
    input_args = {'PlotLayout', false, ... 
                  'PlotPML', false, ...
                  'DisplayMask', source.p_mask | focus_points_mask | model==255,...
                  'RecordMovie', true, 'MovieName', strcat('FocusTarget',num2str(ipoint)),...
                  'DataCast', 'single'};
          
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    if size(target_points,1)>1
        figure;
        p_max = reshape(sensor_data.p_max, Nx, Ny);     
        p_max(model<21 | model>221) = 0;           % ignora pontos fora do cranio
        imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p_max);
        h = colorbar; xlabel(h, '[Pa]');
        title(strcat('Max Acoustic Pressure - Point (', num2str(target_points(ipoint,1)), ', ', num2str(target_points(ipoint,2)), ')'));
    end
    
    % add new sensor data to overal pressure data
    total_sensor_max_pressure = max(total_sensor_max_pressure,sensor_data.p_max);   % max
end


figure;
p_max = reshape(total_sensor_max_pressure, Nx, Ny);     
p_max(model<21 | model>221) = 0;           % ignora pontos fora do cranio
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p_max);
h = colorbar; xlabel(h, '[Pa]');
title('Max Acoustic Pressure');

