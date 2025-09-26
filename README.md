Rocksat summer 2025
Maxx Menzies
Used to determine how much the center of mass will shift as lunasats deploy. If it exceeds the allowed range then it returns a warning.

%adjustable parameters
mass_lunasat = 200; %grams
motor_clearance_diameter = 4; %cm
barrel_diameter = 9; %cm
num_barrels = 3;
num_lunasats_per_barrel = 3;
location_of_com_wo_lunasats = [0,0]; %cm
total_mass = 13607.8; %grams

%determined by the above inputs
radial_com_lunasat = (barrel_diameter/2) + (motor_clearance_diameter/2); %cm
mass_wo_lunasats = total_mass - (num_lunasats_per_barrel * num_barrels * mass_lunasat); %grams
angle_between_lunasats = 360 / num_barrels;
com_wo_lunasats = mass_wo_lunasats .* location_of_com_wo_lunasats;
