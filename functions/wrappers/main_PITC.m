clear;clc;
path2bids='C:\Users\sdem348\Desktop';
params=init_params(path2bids);
params.PltFlag=1;
params.SaveData=1;
Results=struct;
for i=1:20
    subject=strcat('sub-',num2str(i,'%03.0f'));
    params.subject=subject;
    %% wrap below
    fprintf(strcat('Now Processing Case:',num2str(i),'\n'))
    path2data=string(fullfile(params.data_dir,'derivatives\QVT',subject));
    [data_struct,Labels]=import_QVTplusData(path2data);
    [PI_scat]=enc_PITC_process(data_struct,Labels,path2data,params,[]);
    [PITC,globPI] = enc_PITC_fit(PI_scat,path2data,params);
    [DF,PI,LocFlows,LocFlowErr] = enc_PI_DF_Flows(PI_scat,data_struct,Labels,path2data,params);
    [time,Flow,FlowErr]=enc_HQVesselFlows(data_struct,Labels,params);
    plot_MeanVesselFlows(time,Flow,FlowErr)
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