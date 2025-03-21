clear; clc;

% Define paths
% path_to_data = 'E:\PhD\data\CNIC\BMRI198894\';
% output_path = 'C:\Users\u149879\Desktop\autoQVT';
path_to_data = '/media/DATA/afernandezpe/Projects/Reg4DQFlow/CasosNuevos/Data/PESA13593969/';
output_path = '/media/DATA/afernandezpe/Projects/Reg4DQFlow/CasosNuevos/QVTplus/PESA13593969/';

% Load data
[data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

% Perform label transfer and preprocessing
% eICAB_path = 'E:\PhD\data\processed\QVT\BMRI198894';
eICAB_path = '/media/DATA/afernandezpe/Projects/Reg4DQFlow/CasosNuevos/eICAB/PESA13593969/TOF/NIFTI/';
[correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, path_to_data, data_struct);

% Generate LOCs
[correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

% Save vessel-specific data automatically
saveVesselData(LOCs, data_struct, output_path);

%Save data for qvt+
generateQVTplus(correspondenceDict, LOCs, output_path)

disp('Processing completed successfully.');
