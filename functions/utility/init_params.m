function params = init_params(data_dir)
    params=struct;
    params.Version='v0.1.0';
    params.SaveData=1;
    params.PltFlag=1; %Make plots
    params.data_dir=data_dir;
    params.thresh=2.5;
    params.rerun=0;
    params.FlowType='Local';
    params.subject=[];
end