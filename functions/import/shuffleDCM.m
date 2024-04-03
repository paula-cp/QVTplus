%% help function [V,dcminfo] = shuffleDCM(path,flip)
% This function loads and organises the DICOM 4D flow files into a matrix
% in Matlab. 
% % Typically DICOM file organisation is stored in the file name. but may 
% be sorted nonsequentially so this function will find the identifying 
% slice number that tells the order, then stack into the matrix properly. 
%
% flip is in the instance the ordering is reversed, rarely the case, but
% for some reason, happens. Thank you GE! Leave zero unless something is
% wonkey
%
% You may need to add your own regular expression to load the data. 
function [V,dcminfo] = shuffleDCM(path,flip)
DIR=dir(path);
% These are the regular expressions for difference DICOM files I've come
% across, if yours is different, please add it to the cell and the function
% will naturally become more robust. 
Exp={'i\d*.MRDC.(\d*)$';
    '.*-(\d*).dcm$';
    'IMG.(\d*_\d*)';
    '.*i(\d*)$';};
expID=0;
flag=0;
Count=1;
for i=1:length(DIR)
    Name=DIR(i).name;
    if expID==0
        for j=1:length(Exp)
            [match] = regexp(Name,Exp{j},'tokens');
            if length(match) == 1
                expID=j;
                flag=1;
            end
        end
    end
    if flag==1
        [match] = regexp(Name,Exp{expID},'tokens');
        if length(match) == 1
            idx=str2double(match{1});
            Nums(Count,:)=[idx i];
            Count=Count+1;
        end
    end
end
Sorted = sortrows(Nums,1);
for i=1:length(Sorted)
    SortedDir{i,1}=DIR(Sorted(i,2)).name;
end
%At this point, SortedDir gives us the filenames in order to be loaded,
%below just loads the matrix. Stacking is based on the number of timepoints
% (numphase), and slices per stack.

filename=SortedDir{1};
slice = dicomread(fullfile(path,filename));
[a,b]=size(slice);
dcminfo = dicominfo(fullfile(path,filename));
numphase=15;%dcminfo.CardiacNumberOfImages;
slices=length(SortedDir)/numphase;
V=zeros([a,b,slices,numphase],'single');
if flip==0
    pos=0;
else
    pos=slices+1;
end

numphase = 1;
for i=1:length(SortedDir)
    filename=SortedDir{i};
    slice = dicomread(fullfile(path,filename));
    phase = rem(i,(numphase*101)+1);

    if phase ~=0
        pos = pos+1;
        phase = numphase;
    end

    if phase == 0
        pos = 1;
        numphase = numphase + 1;
        phase = numphase;
    end
    [~,name]=fileparts(path);
    if name == 'P3'
        V(:,:,pos,phase)=slice(:,:);
    elseif name == 'P2'
        V(:,:,pos,phase)=slice(:,:);
    else
        V(:,:,pos,phase)=slice(:,:);
    end
    %V(:,:,pos,phase)=slice(:,:);
    dicti{pos,phase}=filename;
end