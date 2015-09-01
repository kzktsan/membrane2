function [allPSFs4D ell2norm2_map]=makeFovPSFs(ell1norm,targetEll2norm2,rho,theta,extraPadPSF,PatchOffset)
%
% Builds the constrained set of blurring PSFs for the foveation operator:         (ver 2.00, Oct. 4, 2014)
% each PSF has ell^1 norm equal to ell1norm and the PFS to be used when blurring
% pixels at coordinates (i,j) has ell^2 norm equal to targetEll2norm2(i,j)
%
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


if ~exist('extraPadPSF','var')||isempty(extraPadPSF);
    extraPadPSF=0;
end
if ~exist('rho','var')||isempty(rho);
    rho=1;
end
if ~exist('theta','var')||isempty(theta);
    theta=0;
end
if ~exist('PatchOffset','var')||isempty(PatchOffset);
    PatchOffset=[0 0];
end

[size_U_1,size_U_2]=size(targetEll2norm2); % size of the ell-2 norm map; coincides with the size of k_kernel and size of the foveated patch.

% upper and lower bounds for the squared l2norm of the blurring PSF (needed to ...... guarantee foveation operator properties in [1]? TBC
% TBD: is this l2norm or l2norm squared?  CHECK ell2norm_Lbound!!
ell2norm_Ubound=ell1norm^2;   % this is attained for a Dirac-delta impulse
ell2norm_Lbound=ell1norm^2/((size_U_2+2*extraPadPSF)*(size_U_1+2*extraPadPSF));  % this is attained for a uniform PSF over the whole PSF support.
targetEll2norm2=min(ell2norm_Ubound-sqrt(eps),max(targetEll2norm2,ell2norm_Lbound+sqrt(eps)));

ell2norm2_map=zeros(size_U_1,size_U_2);
allPSFs4D =zeros(size_U_1+extraPadPSF*2,size_U_2+extraPadPSF*2,size_U_1,size_U_2);

varsigmaInit=ell1norm./(2*sqrt(pi*targetEll2norm2));
if any(extraPadPSF)
    disp(' - Computing blurring PSFs to be used by the foveation operator ...  ')
    textString=['extraPadPSF=',num2str(extraPadPSF)];
else
    disp(' - Computing "patch-supported" blurring PSFs to be used by the foveation operator ...  ')
    textString=[];
end
if any(PatchOffset);
    textString=[textString,'   PatchOffset=[',num2str(PatchOffset(1)),',',num2str(PatchOffset(2)),']'];
end
if any(PatchOffset)||any(extraPadPSF)
    disp(['     ',textString]);
end

timeFirst=now;
timeOld=timeFirst;
loopCount=0;
loopCountPartial=0;
loopTotal=size_U_1*size_U_2;
progressStringEveryNseconds=2.0;
timeLoop=inf;
for u_1=1:size_U_1;
    for u_2=1:size_U_2;
        
        % define cost function fitCostFun  which penalizes the difference between ell2norm2 of the PSF and targetEll2norm2.
        fitCostFun=@(ell2norm2,targetEll2norm2) abs(ell2norm2-targetEll2norm2);
        fTol=targetEll2norm2(u_1,u_2)*eps('single'); % numerical tolerance of fitCostFun to decide convergence
        MaxFunEvals=40; % maximum number of evaluations of the cost function
        showFig=0; % show PSF figure during optimization?
        varsigma=varsigmaInit(u_1,u_2); % value of varsigma used as initialization for the iterative optimization
        
        if 1 % use iterative scaling of varsigma
            fitError=+inf;
            funEvalCount=0;
            while fitError>fTol
                funEvalCount=funEvalCount+1;
                ell2norm2=blurEll2norm2(u_1,u_2,size_U_1,size_U_2,extraPadPSF,varsigma,ell1norm,rho,theta,showFig,PatchOffset);
                fitError=fitCostFun(ell2norm2,targetEll2norm2(u_1,u_2));
                if fitError>fTol
                    ratioNew=sqrt(ell2norm2/targetEll2norm2(u_1,u_2)); % varsigma is (very roughly) inversely proportional to sqrt(ell2norm2), therefore this is good guess of scaling factor for correcting varsigma
                    if funEvalCount>1
                        ratioPow=ratioPow*log(ratioOld)/log(ratioOld/ratioNew);  % here we try to guess how many times one should actually have to multiply varsigma by ratioNew to get the ratio sqrt(ell2norm2/targetEll2norm2(u_1,u_2)) to equal one. ratioPow is approximately 1 whenever the above inverse proportionality holds accurately. ratioPow becomes very large when targetEll2norm2 is close to the its upper limiting value, i.e. when the PSF should become a Dirac.  This provides a substantial acceleration of the convergence.
                    else
                        ratioPow=1;
                    end
                    % disp(num2str(ratioPow));
                    varsigma=varsigma*ratioNew^ratioPow;
                    ratioOld=ratioNew;
                end
                if funEvalCount>=MaxFunEvals
                    break
                end
            end
        else % varsigma optimization using Nelder-Mead algorithm (fminsearch)
            MaxIter=ceil(MaxFunEvals/2);  % maximum number of iterations of fminsearch
            % optimSetDisplay='iter';
            optimSetDisplay='none';
            xTol=varsigma*eps('single');  % numerical tolerance of varsigma to decide convergence
            
            % use fminsearch to find correct value of varsigma (it is necessary to optimize varsigma if extraPadPSF=0, since blur PSFs may have a support that does not fit over the actual domain of the desired foveated patch)
            [varsigma fitError]=fminsearch(@(varsigma) fitCostFun(blurEll2norm2(u_1,u_2,size_U_1,size_U_2,extraPadPSF,varsigma,ell1norm,rho,theta,showFig,PatchOffset),targetEll2norm2(u_1,u_2)),varsigma,optimset('MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'TolX',xTol,'TolFun',fTol,'Display',optimSetDisplay));
            
        end
        showFig=0;
        [ell2norm2 blurPSF]=blurEll2norm2(u_1,u_2,size_U_1,size_U_2,extraPadPSF,varsigma,ell1norm,rho,theta,showFig,PatchOffset);
        ell2norm2_map(u_1,u_2)=ell2norm2;
        allPSFs4D (:,:,u_1,u_2)=blurPSF;
        
        
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

% % this is a figure to verify the accuracy of the optimization of varsigma:
% fitError=abs(ell2norm2_map-targetEll2norm2);  figure, imagesc(log10(fitError),[-13 -3]),colorbar;

end

%%


function [ell2norm2,blurPSF]=blurEll2norm2(u_1,u_2,size_U_1,size_U_2,extraPadPSF,varsigma,ell1norm,rho,theta,showFig,PatchOffset)
% this function builds a PSF with a given ell1norm and varsigma value
% OUTPUTS:
%   ell2norm2 squared \ell-2 norm of the PSF
%   blurPSF is the PSF
%
% INPUTS:
%   u_1,u_2   coordinate of the center of the PSF
%   size_U_1,size_U_2    size of the foveated patch
%   extraPadPSF   extra padding of the patch domain
%   ell1norm   ell1norm of the PSF
%   rho,theta    anisotropy parameters
%   PatchOffset  spatial shift of the operator (must be absent or set to [0,0])
%   showFig    useful to display PSFs during optimization


if varsigma<=0  % minimal degenerate case
    
    blurPSF=zeros(2*extraPadPSF+[size_U_1,size_U_2]);  blurPSF(u_1+extraPadPSF,u_2+extraPadPSF)=1;     % discrete Dirac
    
elseif varsigma==inf; % maximal degenerate case
    
    blurPSF=ones(2*extraPadPSF+[size_U_1,size_U_2]);   % uniform PSF
    
else
    
    % center of the patch
    Ucenter_1=(size_U_1+1)/2;
    Ucenter_2=(size_U_2+1)/2;
    
    if exist('PatchOffset','var');
        u_1=u_1+PatchOffset(1);
        u_2=u_2+PatchOffset(2);
        Ucenter_1=Ucenter_1+PatchOffset(1);
        Ucenter_2=Ucenter_2+PatchOffset(2);
    end
    
    
    % oversampling factor
    overSamplingFactor=5;   % use oversampling to mitigate issues with discretization of Gaussian kernels
    overSamplingFactor=ceil((overSamplingFactor-1)/2)*2+1;   % overSamplingFactor must be an odd positive integer
    
    % create oversampled rotated grid with origin at the impulse position
    [grid_2,grid_1]=meshgrid(linspace(-extraPadPSF+1-1/overSamplingFactor,size_U_2+extraPadPSF+1/overSamplingFactor,(2*extraPadPSF+size_U_2)*overSamplingFactor)-u_2,linspace(-extraPadPSF+1-1/overSamplingFactor,size_U_1+extraPadPSF+1/overSamplingFactor,(2*extraPadPSF+size_U_1)*overSamplingFactor)-u_1);
    theta_prime=atan2((Ucenter_2-u_2),(Ucenter_1-u_1))+theta;   %% theta_prime gives the rotation between the cardinal axes of the psf and the cardinal axes of the patch.
    gridRotTan=grid_2*cos(theta_prime)-grid_1*sin(theta_prime);   % this is tangential if theta=0;
    gridRotRad=grid_2*sin(theta_prime)+grid_1*cos(theta_prime);   % this is radial if theta=0;
    
    
    % define standard deviations sigmaTan and sigmaRad of the bivariate Gaussian blurring PSF with respect to the rotated grid.
    radius2=(Ucenter_2-u_2)^2+(Ucenter_1-u_1)^2;
    if radius2==0
        rho=1;   % enforce the circularly-symmetric PSF at the fixation point
    elseif ~exist('rho','var')
        rho=4.0;   % default to anisotropic with rho=4.0
    end
    sigmaTan=varsigma/sqrt(rho);  % this is tangential if theta=0;
    sigmaRad=varsigma*sqrt(rho);  % this is radial if theta=0;
    
    % build un-normalized elliptical Gaussian blurring PSF on oversampled rotated grid
    blurPSF= exp( -( gridRotTan.^2/(2*sigmaTan^2) + gridRotRad.^2/(2*sigmaRad^2) ) )  /  (2*pi*sigmaTan*sigmaRad);
    % downsample to remove oversampling
    blurPSF=conv2(blurPSF,ones(overSamplingFactor)/overSamplingFactor^2,'same');  % basic antialiasing using uniform kernel
    blurPSF=blurPSF((overSamplingFactor-1)/2+1:overSamplingFactor:end-(overSamplingFactor-1)/2,(overSamplingFactor-1)/2+1:overSamplingFactor:end-(overSamplingFactor-1)/2);  % decimation
end

blurPSF=blurPSF./sum(blurPSF(:))*ell1norm;  % enforce l1-norm of blurPSF  (in older codes, this used to be sqrt(-2*pi/log(1-p))*sqrt(kappa_0))
ell2norm2=sum(sum(blurPSF.*blurPSF));       % compute l2-norm^2 of blurPSF


if exist('showFig','var')&&showFig
    figure(5)
    axisHandle=gca;
    image(blurPSF,'parent',axisHandle,'CDataMapping','scaled');colormap(hot)
    title(axisHandle,['u=[',num2str(u_1),',',num2str(u_2),']  \varsigma=',num2str(varsigma),'  squared \ell-2 norm =',num2str(ell2norm2)]);
    colorbar,axis equal tight
    pause(eps)
end

end
