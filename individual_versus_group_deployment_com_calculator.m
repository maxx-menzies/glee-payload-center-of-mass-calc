%% Maxx Menzies
%RockSat Summer 2025
%For questions contact: maxx.menzies@colorado.edu
%Note: XY Axis Center of Mass --> com
close all; clear; clc;

%% Set Up
%adjustable parameters
mass_lunasat = 200; %grams
motor_clearance_diameter = 4; %cm
barrel_diameter = 9; %cm
num_barrels = 3;
num_lunasats_per_barrel = 3;
location_of_com_wo_lunasats = [0,0]; %cm
total_mass = 13607.8; %grams

%determined by the above
radial_com_lunasat = (barrel_diameter/2) + (motor_clearance_diameter/2); %cm
mass_wo_lunasats = total_mass - (num_lunasats_per_barrel * num_barrels * mass_lunasat); %grams
angle_between_lunasats = 360 / num_barrels;
com_wo_lunasats = mass_wo_lunasats .* location_of_com_wo_lunasats;

%preallocations
lunasat_pos = zeros(num_barrels,2);
lunasat_com = zeros(num_barrels*num_lunasats_per_barrel,2);
group_CoM_cm = zeros((num_barrels * num_lunasats_per_barrel) + 1,2); %cm
individual_CoM_cm = zeros((num_barrels * num_lunasats_per_barrel) + 1,2); %cm
Number_Deployed = zeros(((num_barrels * num_lunasats_per_barrel) + 1), 1);
lunasats_deployed_indiv = 0;
lunasats_deployed_group = 0;

%populate the lunasat position array
for j = 1:num_barrels
    lunasat_pos(j,1) = radial_com_lunasat * cosd((j - 1) * angle_between_lunasats);
    lunasat_pos(j,2) = radial_com_lunasat * sind((j - 1) * angle_between_lunasats);
end

%populate the lunasat com array
k = 1;
for i = 1:num_lunasats_per_barrel:num_barrels*num_lunasats_per_barrel
    for j = 0:(num_lunasats_per_barrel - 1)
        lunasat_com(i+j,:) = mass_lunasat * lunasat_pos(k,:);
    end
    k = k + 1;
end


%% Find com data
total_mass0 = total_mass;
lunasat_com0 = lunasat_com;

%eject by group
for i = 1:((num_barrels * num_lunasats_per_barrel) + 1)
    group_CoM_cm(i,:) = (com_wo_lunasats + sum(lunasat_com)) / total_mass;
    lunasat_com(i,:) = 0;
    total_mass = total_mass - mass_lunasat;
end
total_mass = total_mass0;
%repopulate the lunasat_com array
for i = 1:num_barrels:num_barrels*num_lunasats_per_barrel
    for j = 0:num_barrels - 1
        lunasat_com(i+j,:) = mass_lunasat * lunasat_pos(j+1,:);
    end
end


%eject individually
for i = 1:((num_barrels * num_lunasats_per_barrel) + 1)
    individual_CoM_cm(i,:) = (com_wo_lunasats + sum(lunasat_com)) / total_mass;
    lunasat_com(i,:) = 0;
    total_mass = total_mass - mass_lunasat;
end


%% Data Display
%table
for i = 1:((num_barrels * num_lunasats_per_barrel) + 1)
    Number_Deployed(i) = i - 1;
end
data = table(Number_Deployed, individual_CoM_cm, group_CoM_cm);
disp(data);
%alert
[num_indiv_row, num_indiv_col] = size(individual_CoM_cm);
for i = 1:num_indiv_row
    for j = 1:num_indiv_col
        if (individual_CoM_cm(i,j) >= 1.27) %1.27cm = 0.5in
            disp("****Individual Deployment exceeds CoM regulations****");
        end
        if (group_CoM_cm(i,j) >= 1.27) %1.27cm = 0.5in
            disp("****Group deployment exceeds CoM regulations****");
        end
    end
end
%plot
%indiv_radii = hypot(individual_CoM_cm(:,1), individual_CoM_cm(:,2));
%group_radii = hypot(group_CoM_cm(:,1), group_CoM_cm(:,2));
%report
indiv_x_max = max(abs(individual_CoM_cm(:,1))); %cm
indiv_y_max = max(abs(individual_CoM_cm(:,2))); %cm
group_x_max = max(abs(group_CoM_cm(:,1))); %cm
group_y_max = max(abs(group_CoM_cm(:,2))); %cm
indiv_x_tolerance = 1.27 - indiv_x_max; %cm
indiv_y_tolerance = 1.27 - indiv_y_max; %cm
group_x_tolerance = 1.27 - group_x_max; %cm
group_y_tolerance = 1.27 - group_y_max; %cm
fprintf('The limit is +/- 0.5in or 1.27cm in either direction. \n');
fprintf('The individual deployment max x coordinate is %.4fcm (tolerance: %.4fcm) and the max y is %.4fcm (tolerance: %.4fcm) \n', indiv_x_max, indiv_x_tolerance, indiv_y_max, indiv_y_tolerance);
fprintf('The group deployment max x coordinate is %.4fcm (tolerance: %.4fcm) and the max y is %.4fcm (tolerance: %.4fcm) \n', group_x_max, group_x_tolerance, group_y_max, group_y_tolerance);