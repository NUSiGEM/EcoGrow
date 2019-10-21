function response_surface = Surface_Response(selectioninput1,selectioninput2,tmax,scale);
% close all;
% clear all;
format shortE;
global titles tspan;

titles = {'OD' 'Toxin' 'Antitoxin'};
tspan  = [0:30:tmax];

args = containers.Map();
if (scale == 1)
    args(selectioninput1) = {scale [1 2.5 5 7.5]' 1};
    args(selectioninput2) = {scale [1 2.5 5 7.5]' 1};

elseif (scale == 0)
    switch(selectioninput1)
        case {"mu_max"}
          magscale1 = 0.01;
        case {"OD_max"}
          magscale1 = 0.001;
        case {"OD_min"}
          magscale1 = 0.01;
        case {"nOD"}
          magscale1 = 0.1;
        case {"Ktox"}
          magscale1 = 1e-5;
        case {"ntox"}
          magscale1 = 0.01;
        case {"Kb"}
          magscale1 = 1e-6;
        case {"synt"}
          magscale1 = 1e-7;
        case {"Ki"}
          magscale1 = 1e-3;
        case {"degt"}
          magscale1 = 1e-6;
        case {"syna"}
          magscale1 = 1e-7;
        case {"Ka"}
          magscale1 = 1e-1;     
        case {"dega"}
          magscale1 = 1e-6;
        case {"ti"}
          magscale1 = 1;
        case {"IPTG"}
          magscale1 = 1e-1;
        case {"ta"}
          magscale1  = 1e2;
        case {"Ara"}
          magscale1 =1e-1;
        case {"OD_init"}
          magscale1 = 1e-2;
    end 
    switch(selectioninput2)
        case {"mu_max"}
          magscale2 = 0.01;
        case {"OD_max"}
          magscale2 = 0.001;
        case {"OD_min"}
          magscale2 = 0.01;
        case {"nOD"}
          magscale2 = 0.1;
        case {"Ktox"}
          magscale2 = 1e-5;
        case {"ntox"}
          magscale2 = 0.01;
        case {"Kb"}
          magscale2 = 1e-6;
        case {"synt"}
          magscale2 = 1e-7;
        case {"Ki"}
          magscale2 = 1e-3;
        case {"degt"}
          magscale2 = 1e-6;
        case {"syna"}
          magscale2 = 1e-7;
        case {"Ka"}
          magscale2 = 1e-1;     
        case {"dega"}
          magscale2 = 1e-6;
        case {"ti"}
          magscale2 = 1;
        case {"IPTG"}
          magscale2 = 1e-1;
        case {"ta"}
          magscale2  = 1e2;
        case {"Ara"}
          magscale2 =1e-1;
        case {"OD_init"}
          magscale2 = 1e-2;
    end  
   args(selectioninput1) = {scale magscale1*[-1 0 1 2 3 4 5 6 7 8]' 1};
   args(selectioninput2) = {scale magscale2*[-1 0 1 2 3 4 5 6 7 8]' 1};
end 

switch(selectioninput1)
    case {"mu_max"}
      selectioncol1 = 1; 
    case {"OD_max"}
      selectioncol1 = 2;
    case {"OD_min"}
      selectioncol1 = 3;
    case {"nOD"}
      selectioncol1 = 4; 
    case {"Ktox"}
      selectioncol1 = 5;
    case {"ntox"}
      selectioncol1 = 6;
    case {"Kb"}
      selectioncol1 = 7; 
    case {"synt"}
      selectioncol1 = 8;
    case {"Ki"}
      selectioncol1 = 9;
    case {"degt"}
      selectioncol1 = 10; 
    case {"syna"}
      selectioncol1 = 11;
    case {"Ka"}
      selectioncol1 = 12;     
    case {"dega"}
      selectioncol1 = 13;
    case {"ti"}
      selectioncol1 = 14;
    case {"IPTG"}
      selectioncol1 = 15;
    case {"ta"}
      selectioncol1 = 16; 
    case {"Ara"}
      selectioncol1 = 17;
    case {"OD_init"}
      selectioncol1 = 18;
end  

switch (selectioninput2)
        case {"mu_max"}
      selectioncol2 = 1; 
    case {"OD_max"}
      selectioncol2 = 2;
    case {"OD_min"}
      selectioncol2 = 3;
    case {"nOD"}
      selectioncol2 = 4; 
    case {"Ktox"}
      selectioncol2 = 5;
    case {"ntox"}
      selectioncol2 = 6;
    case {"Kb"}
      selectioncol2 = 7; 
    case {"synt"}
      selectioncol2 = 8;
    case {"Ki"}
      selectioncol2 = 9;
    case {"degt"}
      selectioncol2 = 10; 
    case {"syna"}
      selectioncol2 = 11;
    case {"Ka"}
      selectioncol2 = 12;     
    case {"dega"}
      selectioncol2 = 13;
    case {"ti"}
      selectioncol2 = 14;
    case {"IPTG"}
      selectioncol2 = 15;
    case {"ta"}
      selectioncol2 = 16; 
    case {"Ara"}
      selectioncol2 = 17;
    case {"OD_init"}
      selectioncol2 = 18;
end 
T = readtable('20_runs.csv', 'ReadRowNames', true);
t = T(1,:);

p    = permute_params(t, args, 1);


    response_surface = get_response_surface(p, @model_hic, tspan, [selectioncol1 selectioncol2], 0);


    response_surface = array2table(response_surface, 'VariableNames', {selectioninput1,selectioninput2, 'Objective'});


end

function response_surface = get_response_surface(p, model, tspan, indices, make_plot)
index1 = indices(1);
index2 = indices(2);

response_surface = zeros(length(p), 3);
data             = containers.Map('KeyType','int32','ValueType','any');
if make_plot == 1
    fig      = figure;
    AX{1, 3} = [];
    for i =1:3
        AX{i} = subplot(1, 3, i);
    end
else
    AX = 0;
end

size_p   = size(p);
height_p = size_p(1);
for i=1:height_p
    params = p(i,:);
    init   = [params(18) 0 0];
    y_model = plot_model_1(model, tspan, init, params, AX);
    data(i) = [y_model(:,1,1) y_model(:,2, 1) y_model(:,2, 2) y_model(:,2, 3)]; 

    response_surface(i,:) = [params(index1) params(index2) get_objective(data(i), params)];
end

end

function result = get_objective(dataset, params)
recovery   = get_recovery(dataset, params);
reduction  = get_reduction(dataset, params);
IPTG_usage = get_IPTG_usage(dataset, params);
Ara_usage  = get_Ara_usage(dataset, params);
arrest     = is_arrest(dataset, params); 

result = arrest*(reduction)*recovery*(IPTG_usage^0.5)*(Ara_usage^0.5);%

end

function result = is_arrest(dataset, params)

add_IPTG_after = params(14)+120;
add_Ara        = params(16);

index1 = find(dataset(:,1) >= add_IPTG_after, length(dataset(:,1)), 'first');
index2 = find(dataset(:,1) >= add_Ara, length(dataset(:,1)), 'first');

OD1 = dataset(index1(1),2);
OD2 = dataset(index2(1),2);

result = (OD2 - OD1)/OD1 <= 0.1;

end

function result = get_reduction(dataset, params)

add_IPTG_after = params(14)+120;
add_Ara        = params(16);

index1 = find(dataset(:,1) >= add_IPTG_after, length(dataset(:,1)), 'first');
index2 = find(dataset(:,1) >= add_Ara, length(dataset(:,1)), 'first');

OD1 = dataset(index1(1),2);
OD2 = dataset(index2(1),2);

result = ( (OD2 - OD1)/OD1 )^-1;

end

function recovery = get_recovery(dataset, params)
add_Ara = params(16);
index   = find(dataset(:,1) >= add_Ara, length(dataset(:,1)), 'first');

OD1      = dataset(index(1),2);
final_OD = dataset(length(dataset),2);
recovery = final_OD - OD1;
end

function IPTG_usage = get_IPTG_usage(dataset, params)
IPTG       = params(15);
IPTG_usage = min(1, 1/IPTG);
end

function Ara_usage = get_Ara_usage(dataset, params)
Ara       = params(17);
Ara_usage = min(1, 0.2/Ara);
end

function result = collapse_map(map)
result = [];
keyset = keys(map);

for i = 1:length(keyset)
    key    = keyset{i};
    result = [result; map(key)];
end
end

function y_model = plot_model_1(model, tspan, init, params, AX)

num_state_variables = length(init(1,:));
num_init_conditions = length(init(:,1));
num_time_points     = length(tspan);

y_model = zeros(num_time_points,  num_init_conditions+1, num_state_variables);

for i=1:num_state_variables
    y_model(:,1,i) = tspan;
end

for i=1:num_init_conditions
    [t, y] = ode23s(@(t,y) model(t, y, params), tspan, init(i,:));
    for ii=1:num_state_variables
        y_model(:,i+1,ii) = y(:,ii); 
    end
end


if class(AX) == "cell" && length(AX) > 0
   make_plot(y_model, AX);
end

end

function AX1 = make_plot(data, AX)
global titles
if class(AX) == "cell" && length(AX) == 0
    fig       = figure;
    AX1{1, 3} = [];
    for i =1:3
        AX1{i} = subplot(1, 3, i);
    end
else
    AX1 = AX;
end

for i=1:3
        axes(AX{i});
        plot(data(:,1,i), data(:,2,i),'-','linewidth',2); hold on;
        title(titles{i}, 'FontSize', 20);
            
end
    
end

function data = import_params_csv(files, positions)
data = containers.Map;

for i=1:length(files)
    file       = files(i);
    position   = positions{i};
    if length(position) == 2
        M = dlmread(file, ',' ,position(1), position(2));
    else
        M = dlmread(file, ',' ,position);
    end
    
    data(file) = M;
end

end

function dydt = model_hic(t,y,params)
OD    = y(1);
T     = y(2);
A     = y(3);

mu_max = params(1);
OD_max = params(2);
OD_min = params(3);
nOD    = params(4);
Ktox   = params(5);
ntox   = params(6);
Kb     = params(7);
synt   = params(8);
Ki     = params(9);
degt   = params(10);
syna   = params(11);
Ka     = params(12);
dega   = params(13);
ti     = params(14);
IPTG   = params(15);
ta     = params(16);
Ara    = params(17);

if t < ti
    IPTG = 0;
end

if t < ta
    Ara = 0;
end

Tf     = T*Kb/(A + Kb);
mu     =  mu_max*(1-(OD/OD_max)^nOD)*Ktox^ntox /(Ktox^ntox + Tf^ntox);

dOD    = (OD + OD_min)*mu;
dT     = -degt*T + synt*IPTG/(IPTG + Ki);
dA     = -dega*A + syna*Ara /(Ara  + Ka);

dydt = [dOD; dT; dA];
end
