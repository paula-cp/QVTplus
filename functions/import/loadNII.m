function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean,segmentFullEx,autoFlow,pixelSpace, VoxDims, PIvel_val] = loadNII(directory,handles)
% loadNII loads NIFTI files and then processes them.
%
% It returns too much to discuss, but it basically passes through all the
% processed data to paramMap.
%
% Outputs: Everything?
%
% Used by: autoCollectFlow.m, and any separate functions to compute any
% saved data (PITC codes, Damping codes etc)
%
% Notes: The DICOM data may be rotated 180 (reverse flows), this may be
% allieviated by rotating the matrices using something like 
% %MAG=imrotate(MAG(:,:,:),rotImAngle);
% This should be done for all loaded matrices, HOWEVER, this was only
% required on Siemens data for some reason and may not always be needed.
%
% Dependencies: ShuffleDCM.m and all QVT processing codes
%% Initialization
clc
BGPCdone=0; %0=do backgroun correction, 1=don't do background correction.
%VENC = 800; %may change depending on participant
autoFlow=1; %if you want automatically extracted BC's and flow profiles 0 if not.

filetype = 'nii';
set(handles.TextUpdate,'String','Loading .NII Data'); drawnow;

folders = dir(fullfile(directory,'scans'));
folders(ismember({folders.name}, {'.', '..'})) = [];
folders = {folders([folders.isdir]).name};

folderAP = folders(~cellfun('isempty', regexp(folders, 'AP')));
folderAP = folderAP{1};
folderRL = folders(~cellfun('isempty', regexp(folders, 'RL')));
folderRL = folderRL{1};
folderFH = folders(~cellfun('isempty', regexp(folders, 'FH')));
folderFH = folderFH{1};


json_mag = dir(fullfile(directory,'scans', folderAP,'NIFTI','*.json'));
json_mag = fileread(fullfile(json_mag(1).folder,json_mag(1).name));
json_mag = jsondecode(json_mag); 

magvol = dir(fullfile(directory,'scans', folderAP, 'NIFTI', '*.nii.gz'));
magvol = spm_vol(fullfile(magvol(1).folder,magvol(1).name));
mag = flip(spm_read_vols(magvol),3);

vxvol = dir(fullfile(directory,'scans', folderAP, 'NIFTI', '*_ph.nii.gz'));
vxvol = spm_vol(fullfile(vxvol(1).folder,vxvol(1).name));
vx = flip(spm_read_vols(vxvol),3);

vyvol = dir(fullfile(directory,'scans', folderRL, 'NIFTI', '*_ph.nii.gz'));
vyvol = spm_vol(fullfile(vyvol(1).folder,vyvol(1).name));
vy = flip(spm_read_vols(vyvol),3);

vzvol = dir(fullfile(directory,'scans', folderFH, 'NIFTI', '*_ph.nii.gz'));
vzvol = spm_vol(fullfile(vzvol(1).folder,vzvol(1).name));
vz = flip(spm_read_vols(vzvol),3);

[a,c,b,d] = size(vx);
v = zeros([a,c,b,3,d],'single');

% velocities are in cm/s, convert to mm/s

v(:,:,:,2,:)=-squeeze(vx(:,:,:,:))*10;
v(:,:,:,1,:)=squeeze(vy(:,:,:,:))*10;
v(:,:,:,3,:)=-squeeze(vz(:,:,:,:))*10;

set(handles.TextUpdate,'String','Loaded NIfTI data'); drawnow;

vMean = mean(v,5);

MAG = mean(mag,4);

nframes = length(magvol);

% are these needed?
VENC=700;
timeres = 60/(60*nframes);%INFO.NominalInterval/nframes; %temporal resolution (ms) ????

matrix = magvol(1).dim;
VoxDims = sqrt(sum(magvol(1).mat(1:3,1:3).^2));
res = VoxDims(1);
slicespace = VoxDims(3);

%% Import Complex Difference
set(handles.TextUpdate,'String','Loading Complex Difference Data'); drawnow;
timeMIP = calc_angio(MAG, vMean, VENC);

%% Manual Background Phase Correction (if necessary)
back = zeros(size(vMean),'single');
if ~BGPCdone
    set(handles.TextUpdate,'String','Phase Correction with Polynomial'); drawnow;
    [poly_fitx,poly_fity,poly_fitz] = background_phase_correction(MAG,vMean(:,:,:,1),vMean(:,:,:,2),vMean(:,:,:,3));
    disp('Correcting data with polynomial');
    xrange = single(linspace(-1,1,size(MAG,1)));
    yrange = single(linspace(-1,1,size(MAG,2)));
    zrange = single(linspace(-1,1,size(MAG,3)));
    [Y,X,Z] = meshgrid(yrange,xrange,zrange);
    % Get poly data and correct average velocity for x,y,z dimensions
    back(:,:,:,1) = single(evaluate_poly(X,Y,Z,poly_fitx));
    back(:,:,:,2) = single(evaluate_poly(X,Y,Z,poly_fity));
    back(:,:,:,3) = single(evaluate_poly(X,Y,Z,poly_fitz));
    vMean = vMean - back;
    for f=1:nframes
        v(:,:,:,:,f) = v(:,:,:,:,f) - back;
    end 
    clear X Y Z poly_fitx poly_fity poly_fitz xrange yrange zrange
end
%% Find optimum global threshold for total branch segmentation
set(handles.TextUpdate,'String','Segmenting and creating Tree'); drawnow;
step = 0.001; %step size for sliding threshold
UPthresh = 0.8; %max upper threshold when creating Sval curvature plot
SMf = 10;
shiftHM_flag = 1; %flag to shift max curvature by FWHM
medFilt_flag = 1; %flag for median filtering of CD image
own_seg = 0;
if own_seg == 0
    [~,segment] = slidingThreshold(timeMIP,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
    areaThresh = round(sum(segment(:)).*0.005); %minimum area to keep
    conn = 6; %connectivity (i.e. 6-pt)
    segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes
else
    NII = spm_vol('/home/afernandezpe/borrar/4DFlow_reg/tests_QVT/rBMRI198894_4DQflowNeuro_AP_1803_tAvgProj_eICAB_WB.nii.gz');
    segment = spm_read_vols(NII);
    segment(isnan(segment)) = 0;
end
% save raw (cropped) images to imageData structure (for Visual Tool)
imageData.MAG = MAG;
imageData.CD = timeMIP; 
imageData.V = vMean;
imageData.Segmented = segment;
imageData.Header = json_mag;

%save("C:\Users\u149879\Desktop\trial\original_ABI\ABI_nifti",imageData,'-mat')

%% Feature Extraction
% Get trim and create the centerline data
sortingCriteria = 3; %sorts branches by junctions/intersects 
spurLength = 8; %minimum branch length (removes short spurs)
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

%% You can load another segmentation here if you want which will overlap on images
Exseg=segment; %for now, dummy copy, can do feature extraction
%[~,~,branchList2,~] = feature_extraction(sortingCriteria,spurLength,vMean,logical(JSseg),[]);

%% SEND FOR PROCESSING
% Flow parameter calculation, bulk of code is in paramMap_parameters.m
SEG_TYPE = 'thresh'; %kmeans or thresh
if strcmp(SEG_TYPE,'kmeans')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_kmeans(filetype,branchList,matrix,timeMIP,vMean, ...
    back,BGPCdone,directory,nframes,res,MAG, v,slicespace, handles);
elseif strcmp(SEG_TYPE,'thresh')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes,pixelSpace,...
        segmentFullEx,PIvel_val] ...
        = paramMap_params_threshS(filetype,branchList,matrix,timeMIP,vMean, ...
    back,BGPCdone,directory,nframes,res,MAG,handles, v,slicespace,Exseg);
else
    disp("Incorrect segmentation type selected, please select 'kmeans' or 'thresh'");
end 
set(handles.TextUpdate,'String','All Data Loaded'); drawnow;
return