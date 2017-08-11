function [latProf] = matRad_calcLatProf2(baseData,rad_distancesSq,sigmaIni,res,radDepthPoints)

% gotta change the way i evaluate the intervals, maybe i can call it from
% multidose

distance = sqrt(max(rad_distancesSq))./res(1)-8;
% [x,y] = meshgrid(-ceil(dist):1:ceil(dist));
x = -ceil(distance):1:ceil(distance);
radDistSq = bsxfun(@plus,x'.^2,x.^2);

[uniqueRD,ia,ic] = unique(round(radDepthPoints,0),'rows','stable');
% uniqueRD = radDepthPoints(ia);

interpSigma1 = interp1(baseData.depths + baseData.offset,...
    baseData.sigma1 ,radDepthPoints);
interpSigma2 = interp1(baseData.depths + baseData.offset,...
    baseData.sigma2 ,radDepthPoints);
interpWeight = interp1(baseData.depths + baseData.offset,...
    baseData.weight ,radDepthPoints);

interpSigma1 = interpSigma1(ia);
interpSigma2 = interpSigma2(ia); 
interpWeight = interpWeight(ia); 
% nanIdx = isnan(interpSigma1);
% 
% interpSigma1(nanIdx) = [];
% interpSigma2(nanIdx) = [];
% interpWeight(nanIdx) = [];
% uniqueRD(nanIdx) = [];


for i = 1:length(uniqueRD)
    sigNarr = (interpSigma1(i)./res(1)).^2 + (sigmaIni./res(1)).^2;
    sigBro = (interpSigma2(i)./res(1)).^2 + (sigmaIni./res(1)).^2;
    GauNarr = exp( -radDistSq ./ (2*sigNarr))./(2*pi*sigNarr*res(1)^2);
    GauBro = exp( -radDistSq ./ (2*sigBro))./(2*pi*sigBro*res(1)^2);
    latProf.Mask(:,:,i) = baseData.LatCutOff.CompFac .* ((1-interpWeight(i)).*...
        GauNarr + interpWeight(i).* GauBro);
end
disp('')