clear; clc;

% Define paths
path_to_data = 'E:\PhD\data\CNIC\BMRI198894\';
output_path = 'C:\Users\u149879\Desktop\autoQVT';

% Load data
[data_struct, imageData] = loadPreprocessedData(path_to_data, output_path);

% Perform label transfer and preprocessing
eICAB_path = 'E:\PhD\data\processed\QVT\BMRI198894';
[correspondenceDict, multiQVT] = performLabelTransfer(eICAB_path, output_path, imageData, data_struct);

% Generate LOCs
[correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT);

% Save vessel-specific data automatically
saveVesselData(LOCs, data_struct, output_path);

%Save data for qvt+
generateQVTplus(correspondenceDict, LOCs, output_path)

disp('Processing completed successfully.');
