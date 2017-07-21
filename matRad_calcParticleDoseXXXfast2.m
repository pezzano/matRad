function dij = matRad_calcParticleDoseXXXfast2(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for particles...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% rotate the ct cube agreeing with beam direction
rotcube = imrotate(ct.cube{1},-stf.gantryAngle,'bicubic','crop');
ct.cube{1} = permute(rotcube,[2 1 3]);
% rotcube = imrotate(rotcube,-stf.couchAngle,'bicubic','crop');
% ct.cube{1} = permute(rotcube,[1 3 2]);

% reset angles
stf.gantryAngle = 0;
stf.couchAngle = 0;
pln.gantryAngles = 0;
pln.couchAngles = 0;
stf.sourcePoint = stf.sourcePoint_bev;
[stf.ray.targetPoint] = stf.ray.targetPoint_bev;

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;
dij.numOfScenarios     = 1;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
if calcDoseDirect
    numOfBixelsContainer = 1;
else
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end
doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) ...
        && strcmp(pln.radiationMode,'carbon')
    
    alphaDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,dij.numOfScenarios);
    for i = 1:dij.numOfScenarios
        dij.mAlphaDose{i}    = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
        dij.mSqrtBetaDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
    end
    
elseif isequal(pln.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
    dij.RBE = 1.1;
    fprintf(['matRad: Using a constant RBE of 1.1 \n']);
end

% Only take voxels inside patient.
% V = [cst{:,4}];
% V = unique(vertcat(V{:}));

V = find(ct.cube{1} ~=0);

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
    load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
    error(['Could not find the following machine file: ' fileName ]);
end

if isfield(pln,'calcLET') && pln.calcLET
    if isfield(machine.data,'LET')
        letDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
        % Allocate space for dij.dosexLET sparse matrix
        for i = 1:dij.numOfScenarios
            dij.mLETDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
        end
    else
        warndlg('LET not available in the machine data. LET will not be calculated.');
    end
end

% generates tissue class matrix for biological optimization
if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) ...
        && strcmp(pln.radiationMode,'carbon')
    
    fprintf('matRad: loading biological base data... ');
    vTissueIndex = zeros(size(V,1),1);
    
    %set overlap priorities
    cst  = matRad_setOverlapPriorities(cst);
    
    for i = 1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(vertcat(cst{i,4}{:}),V,'rows');
        % check if base data contains alphaX and betaX
        if   isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
            % check if cst is compatiable
            if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX')
                
                IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                    ismember(machine.data(1).betaX,cst{i,5}.betaX));
                
                % check consitency of biological baseData and cst settings
                if ~isempty(IdxTissue)
                    vTissueIndex(row) = IdxTissue;
                else
                    error('biological base data and cst inconsistent\n');
                end
            else
                vTissueIndex(row) = 1;
                fprintf(['matRad: tissue type of ' cst{i,2} ' was set to 1 \n']);
            end
        else
            error('base data is incomplement - alphaX and/or betaX is missing');
        end
        
    end
    fprintf('done.\n');
    
    % issue warning if biological optimization not possible
elseif sum(strcmp(pln.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && ~strcmp(pln.radiationMode,'carbon') ||...
        ~strcmp(pln.radiationMode,'protons') && strcmp(pln.bioOptimization,'const_RBExD')
    warndlg([pln.bioOptimization ' optimization not possible with ' pln.radiationMode '- physical optimization is carried out instead.']);
    pln.bioOptimization = 'none';
end

% compute SSDs
stf = matRad_computeSSD(stf,ct);

fprintf('matRad: Particle dose calculation...\n');
counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);
    
    bixelsPerBeam = 0;
    
    % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
    yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
    zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
    coordsV  = [xCoordsV yCoordsV zCoordsV];
    
    % Get Rotation Matrix
    % Do not transpose matrix since we usage of row vectors &
    % transformation of the coordinate system need double transpose
    
    % rotation around Z axis (gantry)
    rotMat_system_T = matRad_getRotationMatrix(pln.gantryAngles(i),pln.couchAngles(i));
    
    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV = coordsV*rotMat_system_T;
    
    rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);
    
    % Calcualte radiological depth cube
    lateralCutoffRayTracing = 50;
    fprintf('matRad: calculate radiological depth cube...');
    radDepthV = matRad_rayTracing(stf(i),ct,V,rot_coordsV,lateralCutoffRayTracing);
    radDepthsMat = zeros(ct.cubeDim);
    radDepthsMat(V) = radDepthV{1};
    
    fprintf('done.\n');
    
    % get indices of voxels where ray tracing results are available
    radDepthIx = find(~isnan(radDepthV{1}));
    
    % limit rotated coordinates to positions where ray tracing is availabe
    rot_coordsV = rot_coordsV(radDepthIx,:);
    
    % Determine lateral cutoff
    fprintf('matRad: calculate lateral cutoff...');
    cutOffLevel = 1;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    fprintf('done.\n');
    
    % find central ray
    centralIdx = find(stf.rayPos_bev(:,1)==0 & stf.rayPos_bev(:,2)==0 & stf.rayPos_bev(:,3)==0);
    centralEnergyIdx = find([machine.data.energy] == stf.ray(centralIdx).energy(1));
    
    
    % I should have done this with every energy. Now it is good like this,
    % will improve in the future.
    maxLateralCutoffDoseCalc = max(machine.data(centralEnergyIdx).LatCutOff.CutOff);
    
    [ix,radialDist_sq] = matRad_calcGeoDists_old(rot_coordsV, ...
        stf(i).sourcePoint_bev, ...
        stf(i).ray(centralIdx).targetPoint_bev, ...
        machine.meta.SAD, ...
        radDepthIx, ...
        maxLateralCutoffDoseCalc);
    
    radDepths = radDepthV{1}(ix);
    
    if cutOffLevel >= 1
        currIx = radDepths <= machine.data(centralEnergyIdx).depths(end) + machine.data(centralEnergyIdx).offset;
    elseif cutOffLevel < 1 && cutOffLevel > 0
        % perform rough 2D clipping
        currIx = radDepths <= machine.data(centralEnergyIdx).depths(end) + machine.data(centralEnergyIdx).offset & ...
            radialDist_sq <= max(machine.data(centralEnergyIdx).LatCutOff.CutOff.^2);
        
%         % peform fine 2D clipping
%         if length(machine.data(centralEnergyIdx).LatCutOff.CutOff) > 1
%             currIx(currIx) = matRad_interp1((machine.data(centralEnergyIdx).LatCutOff.depths + machine.data(centralEnergyIdx).offset)',...
%                 (machine.data(centralEnergyIdx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
%         end
    else
        error('cutoff must be a value between 0 and 1')
    end
    % obtain radiation depths matrix
    %tic
    totalBeamDose = zeros(ct.cubeDim);
    totalNumberOfLateralFractions = 150;
    
    for j = 1:size(stf(i).rayPerEnergy,1)
        
        energyIdx = find([machine.data.energy]==stf.rayPerEnergy{j,1});
                                                 
        % Here i use mean for SSD, will find a better way in the future
        sigmaIni = matRad_interp1(machine.data(energyIdx).initFocus.dist(1,:)',machine.data(energyIdx).initFocus.sigma(1,:)',sum([stf.ray.SSD])./length([stf.ray.SSD]));
        
        [finalWeight, X1, sigma_sub, radius, posx, posz, numOfSub] = ...
                    matRad_calcWeights(sigmaIni, 2, 'circle');
        
        
        
        [TracMat] = matRad_multipleRayTracing(radDepthsMat,stf(i).rayPerEnergy(j,:),...
            pln.isoCenter,[ct.resolution.x ct.resolution.y ct.resolution.z]);
        
        [latProf] = matRad_calcLatProf2(machine.data(centralEnergyIdx),...
            radialDist_sq, sigmaIni, [ct.resolution.x ct.resolution.y ct.resolution.z],...
            totalNumberOfLateralFractions);
        
        if j==18
            disp(':)')
        end
        
        tempBeamDose = totalBeamDose + matRad_calcMultiDose(TracMat,...
            machine.data(energyIdx),latProf,V,totalNumberOfLateralFractions);
        
        % Display progress and update text 
        matRad_progress(j,size(stf(i).rayPerEnergy,1));
        
        % update waitbar 
        waitbar(j/size(stf(i).rayPerEnergy,1),figureWait);
    end
end

totalBeamDose = zeros(ct.cubeDim);
totalBeamDose(V) = tempBeamDose(V);
dij.physicalDose{1} = totalBeamDose;


try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end
