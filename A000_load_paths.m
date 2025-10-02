%% This script load paths required by GMT-OPT

%% Only change the master directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentDir = pwd;
[parentDir, folderName] = fileparts(currentDir);

root_dir = sprintf('%s/',parentDir);
program_dir = sprintf('%s%s/',root_dir,folderName);

meteo_root_dir = sprintf('%sInputs/fluxnet2015/diel_monthly/',program_dir); % flux tower data
site_info_path = sprintf('%sInputs/site_info_with_LAI.mat',program_dir);
CO2_path = sprintf('%sInputs/monthly_CO2.mat',program_dir);
cali_root=sprintf('%sCalibration_result/',program_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add paths for script search.    Do not mannually change here. Only change if necessary.
addpath(genpath(program_dir))