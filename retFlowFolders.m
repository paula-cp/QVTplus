function [Anatpath,APpath,LRpath,SIpath] = retFlowFolders(path2flow,Vendor,varargin)
% retFlowFolders researches for the variable folder naming of 4Dflow data
%
% path2flow: Path to source directory with all flow folders
% Vendor: Input if naming should follow GE convention, or Siemens.
% varargin: You can optionally input the resolution ("0.5" or ".5") if
% there is more than one 4Dflow data set in your sourcedirectory
%
% NOTE: If your folder namings are dramatically different thant the regular
% expression orded FLOW, then Anat, AP, SI, LR, either change the
% expression or rename your folders.
%
% Outputs: [Anatpath,APpath,LRpath,SIpath]
%   paths: Each path returned is the full path to the 4D flow data for each
%   collection direction and magnitude image.
%
% Used by: loadDICOM.m
%
% Dependencies: None    
    TestDir=dir(path2flow);
    RegExp={'.*FLOW.*(\d{1}(\.)?\d{1}).*Anat.*',...
        '.*FLOW.*(\d{1}(\.)?\d{1}).*AP.*',...
        '.*FLOW.*(\d{1}(\.)?\d{1}).*SI.*',...
        '.*FLOW.*(\d{1}(\.)?\d{1}).*LR.*'};
    if ~isempty(varargin)
        res=varargin(1);
        res=res{1};
        splt=strsplit(res,'.');
        if length(splt) == 2
            strcat(splt{1},splt{2});
        end
    else 
        res=[];
        splt=[];
    end
    Anatpath=[];APpath=[];LRpath=[];SIpath=[];

    %% GE Methods
    if strcmp('GE',Vendor)
        for i=1:length(TestDir)
            Exp = TestDir(i).name;
            [match] = regexp(Exp,RegExp{1},'tokens');
            if length(match) == 1
                if strcmp(match{1},res) | strcmp(match{1},splt) 
                    Anatpath=fullfile(path2flow,Exp);
                elseif isempty(Anatpath)
                    Anatpath=fullfile(path2flow,Exp);
                end
            end
            [match] = regexp(Exp,RegExp{2},'tokens');
            if length(match) == 1
                if strcmp(match{1},res) | strcmp(match{1},splt) 
                    APpath=fullfile(path2flow,Exp);
                elseif isempty(APpath)
                    APpath=fullfile(path2flow,Exp);
                end
            end
            [match] = regexp(Exp,RegExp{3},'tokens');
            if length(match) == 1
                if strcmp(match{1},res) | strcmp(match{1},splt) 
                    SIpath=fullfile(path2flow,Exp);
                elseif isempty(SIpath)
                    SIpath=fullfile(path2flow,Exp);
                end
            end
            [match] = regexp(Exp,RegExp{4},'tokens');
            if length(match) == 1
                if strcmp(match{1},res) | strcmp(match{1},splt) 
                    LRpath=fullfile(path2flow,Exp);
                elseif isempty(LRpath)
                    LRpath=fullfile(path2flow,Exp);
                end
            end
        end
    elseif strcmp('Siemens',Vendor)
        Anatpath=fullfile(path2flow,'M1');
        APpath=fullfile(path2flow,'P3');
        LRpath=fullfile(path2flow,'P2');
        SIpath=fullfile(path2flow,'P1');%P3 -> AP, P2 -> RL, P1 -> SI
    end
end