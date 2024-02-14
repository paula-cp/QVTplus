Path2TOF='C:\Users\sdem348\Desktop\3D TOF FS HyperSense_3';
DIR=dir(Path2Mip);
filename=DIR(3).name;
slice = dicomread(fullfile(Path2Mip,filename));
MIP=zeros([size(slice) 2]);
for i=3:length(DIR)
    filename=DIR(i).name;
    slice = dicomread(fullfile(Path2Mip,filename));
    MIP(:,:,2)=slice;
    MIP(:,:,1)=max(MIP,[],3);
end
imagesc(MIP(:,:,1))
colormap(gray)