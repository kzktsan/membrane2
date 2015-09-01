function AllPSFsTiled=makeAllPSFsTiled(allPSFs4D,gapWidth,gapColor,normalizationType)
%                                                                          (ver 2.00, Oct. 4, 2014)
% generate the figure with the mosaic of blurring PSFs that compose the fovation operator.
%
%  REFERENCES:
%  [1] A. Foi and G. Boracchi, "Foveated self-similarity in nonlocal image filtering",
%  Proc. IS&T/SPIE EI2012 - Human Vision and Electronic Imaging XVII, 8291-32,
%  Burlingame (CA), USA, Jan. 2012.
%  [2] A. Foi and G. Boracchi, "Anisotropically Foveated Nonlocal Image Denoising",
%  Proc. IEEE Int. Conf. Image Process. (ICIP 2013), pp. 464-468, Melbourne, Australia, Sep. 2013.
%  [3] A. Foi and G. Boracchi, "Foveated nonlocal self-similarity", preprint, Oct. 2014.
%
%  Author:  Alessandro Foi, Tampere University of Technology, Finland
%           All rights reserved.
% --------------------------------------------------------------------------------------
%

if ~exist('gapWidth','var')
    gapWidth=0;
end
if ~exist('gapColor','var')
    gapColor=0.8;
end

if ~exist('normalizationType','var')
    normalizationType=0;
end

size_U_1=size(allPSFs4D,3);
size_U_2=size(allPSFs4D,4);
size_BlurPSF_1=size(allPSFs4D,1); % this is typically size_U_1+2*extraPadPSF;
size_BlurPSF_2=size(allPSFs4D,2); % this is typically size_U_2+2*extraPadPSF;

AllPSFsTiled=zeros( size_U_1*(size_BlurPSF_1+gapWidth)-gapWidth,size_U_2*(size_BlurPSF_2+gapWidth)-gapWidth)+gapColor;

for jj1=1:size_U_1;
    for jj2=1:size_U_2;
        BlurPSF=squeeze(allPSFs4D(:,:,jj1,jj2));
        if normalizationType==0
            BlurPSF=BlurPSF./max(abs(BlurPSF(:)));
        elseif normalizationType==1
            BlurPSF=BlurPSF./max(abs(BlurPSF(:)))/2+0.5;
        elseif normalizationType==2
            BlurPSF=BlurPSF./max(abs(allPSFs4D(:)))/2+0.5;
        end
        
        AllPSFsTiled((size_BlurPSF_1+gapWidth)*(jj1-1)+1:(size_BlurPSF_1+gapWidth)*(jj1-1)+size_BlurPSF_1,(size_BlurPSF_2+gapWidth)*(jj2-1)+1:(size_BlurPSF_2+gapWidth)*(jj2-1)+size_BlurPSF_2)=BlurPSF;
    end
end
figHandle=figure;
axisHandle=axes('parent',figHandle);
image(AllPSFsTiled,'parent',axisHandle,'CDataMapping','scaled'), title('Blurring PSFs mosaic'),colormap(figHandle,gray);
axis(axisHandle,'equal','tight','off');
truesize(figHandle);
set(figHandle,'color',max(0,min(1,interp1(linspace(min(AllPSFsTiled(:)),max(AllPSFsTiled(:)),size(get(figHandle,'Colormap'),1)),get(figHandle,'Colormap'),gapColor,'linear','extrap'))));
pause(0.001);
