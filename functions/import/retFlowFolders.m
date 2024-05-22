function [Anatpath,APpath,LRpath,SIpath] = retFlowFolders(path2flow,Vendor,varargin)
% retFlowFolders looks for the variable folder naming of 4Dflow data
%
% path2flow: Path to source directory with all DICOM folders for a
% participant
% Vendor: Input if naming should follow GE convention, or Siemens. This is
% under construction until I can get other vendors!
% varargin: You can optionally input the resolution ("0.5" or ".5") if
% there is more than one 4Dflow data set in your sourcedirectory which is
% tagged by resolution
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
    if ~isempty(varargin)
        res=varargin(1);
        res=res{1};
    else 
        res='';
    end
    RegExp={strcat('.*FLOW.*(',res,').*Anat.*'),...
        strcat('.*FLOW.*(',res,').*AP.*'),...
        strcat('.*FLOW.*(',res,').*SI.*'),...
        strcat('.*FLOW.*(',res,').*LR.*')};
    RegExp={'.*cs.*(\d{1}(\.)?\d{1}).*MC.*',...
        '.*flow.*(\d{1}(\.)?\d{1}).*AP.*',...
        '.*flow.*(\d{1}(\.)?\d{1}).*FH.*',...
        '.*flow.*(\d{1}(\.)?\d{1}).*RL.*'};
    RegExp={'.*flow.*AP', ...
        '.*flow.*FH',...
        '.*flow.*RL'};
    Anatpath=[];APpath=[];LRpath=[];SIpath=[];
    %% GE Methods
    if strcmp('GE',Vendor) | strcmp('Philips',Vendor)
        for i=1:length(TestDir)
            Exp = TestDir(i).name;

            [match] = regexp(Exp,RegExp{1},'tokens');
            if length(match) == 1
                if strcmp(match{1},res) 
                    Anatpath=fullfile(path2flow,Exp,'DICOM');
                elseif isempty(Anatpath)
                    Anatpath=fullfile(path2flow,Exp,'DICOM');
                end
            end
            [match] = regexp(Exp,RegExp{1},'tokens');
            if length(match) == 1
                if strcmp(match{1},res)
                    APpath=fullfile(path2flow,Exp,'DICOM');
                elseif isempty(APpath)
                    APpath=fullfile(path2flow,Exp,'DICOM');
                end
            end
            [match] = regexp(Exp,RegExp{2},'tokens');
            if length(match) == 1
                if strcmp(match{1},res)
                    SIpath=fullfile(path2flow,Exp,'DICOM');
                elseif isempty(SIpath)
                    SIpath=fullfile(path2flow,Exp,'DICOM');
                end
            end
            [match] = regexp(Exp,RegExp{3},'tokens');
            if length(match) == 1
                if strcmp(match{1},res)
                    LRpath=fullfile(path2flow,Exp,'DICOM');
                elseif isempty(LRpath)
                    LRpath=fullfile(path2flow,Exp,'DICOM');
                end
            end
        end
    elseif strcmp('Siemens',Vendor)
        Anatpath=fullfile(path2flow,'M1');
        APpath=fullfile(path2flow,'P3');
        LRpath=fullfile(path2flow,'P2');
        SIpath=fullfile(path2flow,'P1');%P3 -> AP, P2 -> RL, P1 -> SI
    elseif strcmp('Philips',Vendor)

    end
end