function [y_hat]=AnisFovNLM(z,sigma,search_radius,U_radius,disableFov,selfMap,rho,theta,extraPadPSF)
%
%
%  Foveated Nonlocal Means denoising algorithm (ver 2.01, Oct. 9, 2014)
%  Author:  Alessandro Foi, Tampere University of Technology, Finland
% -------------------------------------------------------------------------
%  SYNTAX
%
%     [y_hat] = AnisFovNLM(z,sigma,search_radius,U_radius,disableFov,selfMap,rho,theta,extraPadPSF)
%
%  y_hat :          output denoised image
%
%  z :              input noisy image to be filtered   [REQUIRED]
%  sigma :          standard deviation of the additive white Gaussian noise  [REQUIRED]
%
%  Optional input parameters ---------------------------------------------------------
%  search_radius :  radius of search window     (e.g., 5 <= search_radius <= 13 )
%  U_radius :       radius of patch             (e.g., 1 <= U_radius <= 10)
%  disableFov :     disables foveation when set to 1  (default 0)
%  selfMap :        enables self-map foveation operators when set to 1  (default 1)
%  rho :            elongation of the blur PSFs (rho=1 gives circular-symmetric PSFs)
%  theta :          orientation of the blur PSFs relative to the meridian lines
%  extraPadPSF :    extra padding of the PSFs support
%
%
%
%  This software implements the Foveated Nonlocal Means algorithm
%  proposed in the papers:   (all available online at http://www.cs.tut.fi/~foi/)
%  [1] A. Foi and G. Boracchi, "Foveated self-similarity in nonlocal image filtering",
%  Proc. IS&T/SPIE EI2012 - Human Vision and Electronic Imaging XVII, 8291-32,
%  Burlingame (CA), USA, Jan. 2012.
%  [2] A. Foi and G. Boracchi, "Anisotropically Foveated Nonlocal Image Denoising",
%  Proc. IEEE Int. Conf. Image Process. (ICIP 2013), pp. 464-468, Melbourne, Australia, Sep. 2013.
%  [3] A. Foi and G. Boracchi, "Foveated nonlocal self-similarity", preprint, Oct. 2014.
%
%
%  By setting disableFov=1, this software can replicate the results of the
%  NLmeansfilter.m  code by J.V. Manjon Herrera & A. Buades (31.07.2008),
%  which implements the NL-means filter proposed by A. Buades, B. Coll and
%  J.M. Morel in the paper "A non-local algorithm for image denoising".
%  http://www.mathworks.com/matlabcentral/fileexchange/13176-non-local-means-filter
%
%%%


%% get and verify image size
[size_z_1 size_z_2 numColorChan]=size(z);
if numColorChan>1
    error('Error: current code works for grayscale images only.')
end


%% initialize some variables

% toggles: fovation/windowing, self-map foveation operator
if ~exist('disableFov','var')||isempty(disableFov)
    disableFov=0;  % use foveation by default
end
if ~exist('selfMap','var')||isempty(selfMap)
    selfMap=1;  % use self-map operators by default
end

% anisotropy
if ~exist('theta','var')||isempty(theta)
    theta=0.0*pi;  % theta is the angle between the first axis of the elliptical pdf (main axis if rho>1) and the radial (meridian) line between the psf center and the center of the patch
end
if ~exist('rho','var')||isempty(rho)
    rho=4.0;   % if theta=0, this is the tangent-radial ratio: set to 1 for circular-symmetric blur PSFs, >1 for radially elongated elliptical PSFs, <1 for tangentially elongated elliptical PSFs
end

% By leaving search_radius and U_radius empty, their values are chosen automatically based on the value of sigma (the signal is assumed to be in the standard [0,255] range)
% search-neighborhood and patch radiuses
if ~exist('search_radius','var')||isempty(search_radius)
    sigmaS=[10 20 30 50 70];
    search_radiusS=[8 8 8 8 8]*~disableFov+[5 4 4 5 6]*~~disableFov;
    search_radius=round(interp1(sigmaS,search_radiusS,sigma,'linear','extrap'));
end
if ~exist('U_radius','var')||isempty(U_radius)
    sigmaS=[10 20 30 50 70];
    U_radiusS=[3 5 6 8 9]*~disableFov+[2 5 5 6 7]*~~disableFov;
    U_radius=round(interp1(sigmaS,U_radiusS,sigma,'linear','extrap'));
end


disp(' ');
if ~disableFov
    disp('Computing Foveated NL-means estimate (using foveated patch distance) ...');
else
    disp('Computing NL-means estimate (using windowed patch distance) ...');
end

disp([' . sigma=',num2str(sigma), ', image=',num2str(size_z_2),'x',num2str(size_z_1),', U_radius=',num2str(U_radius),', search_radius=',num2str(search_radius)]);
if disableFov
    textString='windowing (foveation is disabled)';
else
    if rho~=1
        textString=['anisotropic foveation: rho=',num2str(rho), ', theta=',num2str(theta)];
    else
        textString='isotropic foveation';
    end
    if selfMap
        textString=[textString,'  | self-map'];
    end
end
disp([' . ',textString]);


%% Windowing k_kernel  (within the Foveated NL-means this kernel is not used for windowing, but to determine the standard-deviation of the blurring kernels used by the foveation operator)
[u2,u1]=meshgrid(-U_radius:U_radius,-U_radius:U_radius);
if 1 %% kernel used in the NL-means implementation by J.V. Manjon Herrera & A. Buades
    ellInftyDistance=max(abs(u1),abs(u2));
    K_tilde=cumsum((2*(U_radius:-1:0)+1).^-2);
    K_tilde=K_tilde([U_radius,U_radius:-1:1]);
    k_kernel=interp1(0:U_radius,K_tilde,ellInftyDistance,'linear','extrap');
elseif 1 %% a-la Wertheim  (this is just to give an illustration of a different design of the window)
    l2distance=sqrt(u1.^2+u2.^2);
    k_kernel=1-(l2distance*(0.98/U_radius/sqrt(2))).^0.3;
else
    kaiserBeta=8;
    k_kernel=kaiser(U_radius*2+1,kaiserBeta)*kaiser(U_radius*2+1,kaiserBeta)';
    k_kernel=k_kernel.^0.5;
end
k_kernel=k_kernel/sum(k_kernel(:));

numel_k_kernel=numel(k_kernel);



%% Patching
if disableFov
    %% use windowing instead of foveation
    FovMat=diag(sqrt(k_kernel(:)));  % the foveation operator equivalent to windowing is composed by scaled Dirac impulses, hence in matrix form it is diagonal
    if ~exist('extraPadPSF','var')||isempty(extraPadPSF)
        extraPadPSF=0;  % there is no use in having extraPadPSF>0 if foveation is replaced by simple windowing, since the blur PSFs are all Diracs.
    end
else
    %% use foveation instead of windowing
    p=1-exp(-2*pi);
    ell1norm=sqrt(-2*pi/log(1-p))*sqrt(max(k_kernel(:)));  % THIS IS A NUMBER, THE ELL^1 NORM OF ALL THE PSFS, SEE (13) IN [1]
    ell2norm2 = k_kernel; % THIS IS A MATRIX OF THE SIZE OF k_kernel  SEE (14) IN [1]
    
    if selfMap
        extraPadPSF=0;    % in case of self-map foveation operator there is no extra padding of the patch domain.
    else
        if ~exist('extraPadPSF','var')||isempty(extraPadPSF)
            extraPadPSF=ceil(3*sqrt(max(rho,1/rho))*ell1norm./(2*sqrt(pi*min(ell2norm2(:)))));    % use extraPadPSF>0 in order to give more room for the support of the blur PSFs centered near the boundary of the patch (this is a 3sigma rule using the best guess of the maximum standard deviation varsigma of the PSF as ell1norm./(2*sqrt(pi*min(ell2norm2(:))))
        end
    end
    PatchOffset=[0,0];
    [allPSFs4D ~]=makeFovPSFs(ell1norm,ell2norm2,rho,theta,extraPadPSF,PatchOffset); %% construction of blurring PSFs used by the foveation operator with a given ell1norm and ell2norm2
    FovMat=reshape(allPSFs4D,size(allPSFs4D,1)*size(allPSFs4D,2),size(allPSFs4D,3)*size(allPSFs4D,4))';  %% foveation operator is written in matrix form
    % allPSFs4D=reshape(FovMat',U_radius*2+1+extraPadPSF*2,U_radius*2+1+extraPadPSF*2,U_radius*2+1,U_radius*2+1);  %%  reverse operation
    makeAllPSFsTiled(allPSFs4D,2,0.8,0); % generate the figure with the mosaic of blurring PSFs
end

% UcenterIdx=ceil((U_radius*2+1)^2/2);   % index-coordinate of the center of the patch

z_pad=padarray(z,[U_radius+extraPadPSF U_radius+extraPadPSF 0],'symmetric');  % padding (to handle image boundaries without caring about size of patch and blur PSFs support)
BIGfov=zeros(numel_k_kernel,size_z_1,size_z_2); % big 3-D array with all foveated patches;
if disableFov
    for x1=1:size_z_1
        for x2=1:size_z_2
            z_x=z_pad(x1:x1+2*(U_radius+extraPadPSF),x2:x2+2*(U_radius+extraPadPSF));  % extract from the padded input image a patch centered at (x1,x2)
            BIGfov(:,x1,x2)=sqrt(k_kernel(:)).*z_x(:);
        end
    end
else
    for x1=1:size_z_1
        for x2=1:size_z_2
            z_x=z_pad(x1:x1+2*(U_radius+extraPadPSF),x2:x2+2*(U_radius+extraPadPSF));  % extract from the padded input image a patch centered at (x1,x2)
            BIGfov(:,x1,x2)=FovMat*z_x(:);
        end
    end
end

if ~disableFov
    disp(' - Filtering with Foveated NL-means ...  ')
else
    disp(' - Filtering with NL-means (windowing) ...  ')
end

%% Foveated NL-means.

y_hat=zeros(size_z_1,size_z_2);
BIGfov=permute(BIGfov,[2 3 1]);

timeFirst=now;
timeOld=timeFirst;
loopCount=0;
loopCountPartial=0;
loopTotal=size_z_1*size_z_2;
progressStringEveryNseconds=5.0;
timeLoop=inf;
for x1=1:size_z_1
    for x2=1:size_z_2
        
        % set boundaries of search window
        search1min=max(x1-search_radius,1);  search1max=min(x1+search_radius,size_z_1);  % vertical (left; right)
        search2min=max(x2-search_radius,1);  search2max=min(x2+search_radius,size_z_2);  % horizontal (top; bottom)
        search1center=x1-search1min+1; search2center=x2-search2min+1;  % relative coordinates of [x1,x2]
        
        % compute foveated differences and foveated distance for all foveated patches within the search window
        dFov=0;
        for u=1:numel_k_kernel
            dFov=dFov+(BIGfov(search1min:search1max,search2min:search2max,u)-BIGfov(x1,x2,u)).^2; % compute foveated distance
        end
        
        % visualize foveated distances
        if 0
            if ~exist('axisHandle','var')
                figHandle=figure;
                axisHandle=axes('parent',figHandle);
            end
            image((search2min:search2max)-x2,(search1min:search1max)-x1,dFov,'parent',axisHandle,'CDataMapping','scaled');  % colorbar;
            title(axisHandle,num2str([x1 x2]));
            axis(axisHandle,'equal','tight');
            pause(eps)
        end
        
        w=exp(-dFov/(sigma*sigma));
        w(search1center,search2center,:)=0; wmax=max(w(:)); w(search1center,search2center,:)=wmax+~wmax;  % weight for z(x1,x2) is set equal to the maximum weight (unless all weights are zero, in which case it is set to 1)
        w=w/sum(w(:));
        
        % weighted averaging
        y_hat(x1,x2)=sum(sum(w.*z(search1min:search1max,search2min:search2max)));
        
        
        % progress update string
        loopCount=loopCount+1;
        loopCountPartial=loopCountPartial+1;
        if timeLoop*loopCountPartial>progressStringEveryNseconds;
            timeNow=now;
            timeLoop=86400*(timeNow-timeOld)/loopCountPartial;
            if 86400*(timeNow-timeOld)>progressStringEveryNseconds;
                timeOld=timeNow;
                loopCountPartial=0;
                disp(['      ',num2str(floor(100*loopCount/loopTotal)),' %  (ETA ',num2str(timeLoop*(loopTotal-loopCount)),'s)'])
            end
        end
    end
end
disp(['     ',num2str(100),' %  completed  (',num2str(86400*(now-timeFirst)),'s)'])
%% end of code

