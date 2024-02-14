%% Running a Directory
%==================================================
% Running a BIDS processed derivates folder for PI, PI_DF, PI_TC, PWV is below
%==================================================
% This function is designed to run assuming you store your data in BIDS 
% format. If not, modify to run a folder or what not, BIDS is much easier
clear;clc;
path2bids='C:\Users\sdem348\Desktop';
params=init_params(path2bids);
params.PltFlag=1;
params.SaveData=0;
params.rerun=1;
Results=struct;
for i=1:1
    subject=strcat('sub-',num2str(i,'%03.0f'));
    params.subject=subject;
    fprintf(strcat('Now Processing Case:',num2str(i),'\n'))
    path2data=string(fullfile(params.data_dir,'derivatives\QVT',subject));
    [data_struct,Labels]=import_QVTplusData(path2data);
    [PI_scat]=enc_PITC_process(data_struct,Labels,path2data,params,[]);
    [PITC,globPI] = enc_PITC_fit(PI_scat,path2data,params);
    [DF,PI,LocFlows,LocFlowErr] = enc_PI_DF_Flows(PI_scat,data_struct,Labels,path2data,params);
    [time,Flow,FlowErr]=enc_HQVesselFlows(data_struct,Labels,params);
    dirPITC(i,:)=PITC(1,:);
    dirPI(i,:)=[globPI PITC(end,:) PI];
    dirDF(i,:)=DF;
    dirFlow(i,:)=[mean(LocFlows)];
    close all
end
Results.PITC=dirPITC;
Results.PI=dirPI;
Results.DF=dirDF;
Results.Flow=dirFlow;
save(fullfile(params.data_dir,'derivatives\QVT\population\ResultsALL.mat'),"Results")

%% Running 
%==================================================
% Running a single instance of the post processing for PI, PI_DF, PI_TC, PWV is below
%==================================================
clear;clc;
path2data='C:\Users\sdem348\Documents\MATLAB\CHM\QVTplus\testdata'; %should have the  LabelsQVT and qvtData inside
params=init_params(path2data);
params.PltFlag=1; %Plot your results
params.SaveData=1; %Save results and plots
params.rerun=1; %rerun PITCH algorithm, or use stored data
Results=struct;
[data_struct,Labels]=import_QVTplusData(path2data);
[PI_scat]=enc_PITC_process(data_struct,Labels,path2data,params,[]);
[PITC,globPI] = enc_PITC_fit(PI_scat,path2data,params);
[DF,PI,LocFlows,LocFlowErr] = enc_PI_DF_Flows(PI_scat,data_struct,Labels,path2data,params);
[time,Flow,FlowErr]=enc_HQVesselFlows(data_struct,Labels,params);
close all
Results.PITC=PITC(1,:);
Results.PI=[globPI PITC(end,:) PI];
Results.DF=DF;
Results.Flow=[mean(LocFlows)];
save(fullfile(path2data,'ResultsALL.mat'),"Results")
