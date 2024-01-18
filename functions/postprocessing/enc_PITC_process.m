function [PI_scat] = dev_enc_PITC_process(data_struct,Labels,path2data,params)
    %% Load Necessities from labels and processed data
    if ~isempty(path2data)
        if params.rerun == 0 %Don't rerun unless you want to for some reason
            if isfile(fullfile(path2data,'RawPITC.mat'))
                ProcessFlag=0;
                load(fullfile(path2data,'RawPITC.mat'));
            end
        else
            ProcessFlag=1;
        end
    else
        ProcessFlag=1;
    end
    if ProcessFlag==1
        ICA_l = Labels{1,2} ;ICA_l = str2num(ICA_l);
        ICA_r = Labels{2,2} ;ICA_r = str2num(ICA_r);
        BA = Labels{9,2}    ;BA = str2num(BA);
        Exclude=Labels{10,2};Exclude = str2num(Exclude);
        LR=Labels{10,3}     ;LR = str2num(LR);
        startroots = {ICA_l,ICA_r,BA};
        ds = DataStorage();  
        my_app_handle=PITCinteractive(data_struct, startroots, LR, Exclude,ds);
        uiwait(my_app_handle.UIFigure);
        SearchDist=ds.dataArea.SearchDist;
        PI_scat=ds.dataArea.PI_scat;
        G=ds.dataArea.G;
        Mats=ds.dataArea.Mats;
        if params.SaveData==1
            save(fullfile(path2labels,'RawPITC.mat'),'PI_scat','G','Mats','SearchDist');
        end
    end
end