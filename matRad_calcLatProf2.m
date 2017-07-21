function [latProf] = matRad_calcLatProf2(baseData,rad_distancesSq,sigmaIni,res,tot)



distance = sqrt(max(rad_distancesSq))./res(1);
% [x,y] = meshgrid(-ceil(dist):1:ceil(dist));
x = -ceil(distance)-3:1:ceil(distance)+3;
radDistSq = bsxfun(@plus,x'.^2,x.^2);

for i = 1:tot
    ind = baseData.Z <= max(baseData.Z)*i/tot & baseData.Z > max(baseData.Z)*(i-1)/tot;
    if sum(ind)~=0
        if i>1
            bin(i) = floor(mean([find(diff(ind)==1,1)+1, find(diff(ind)==-1,1)]));
            latProf.flag(i) = 1;
        else
            bin(1) = floor(mean([1, find(diff(ind)==-1,1)]));
            latProf.flag(1) = 1;
        end
        sigNarr = (baseData.sigma1(bin(i))./res(1)).^2 + (sigmaIni./res(1)).^2;
        sigBro = (baseData.sigma2(bin(i))./res(1)).^2 + (sigmaIni./res(1)).^2;
%         latProf.Mask(:,:,i) = baseData.LatCutOff.CompFac .* ((1-baseData.weight(bin(i))).*...
%             (exp( -x.^2 ./ (2*sigNarr))./sqrt(2*pi*sigNarr))'*(exp( -x.^2 ./ (2*sigNarr))./sqrt(2*pi*sigNarr)) +...
%             (baseData.weight(bin(i))).*...
%             (exp( -x.^2 ./ (2*sigBro))./sqrt(2*pi*sigBro))'*(exp( -x.^2 ./ (2*sigBro))./sqrt(2*pi*sigBro)));
%         GauNarr = bsxfun(@plus, (exp( -x.^2 ./ (2*sigNarr))./(2*pi*sigNarr))', exp( -x.^2 ./ (2*sigNarr))./(2*pi*sigNarr));
%         GauBro = bsxfun(@plus, (exp( -x.^2 ./ (2*sigBro))./(2*pi*sigBro))', exp( -x.^2 ./ (2*sigBro))./(2*pi*sigBro));
        GauNarr = exp( -radDistSq ./ (2*sigNarr))./(2*pi*sigNarr*res(1)^2);
        GauBro = exp( -radDistSq ./ (2*sigBro))./(2*pi*sigBro*res(1)^2);
        latProf.Mask(:,:,i) = baseData.LatCutOff.CompFac .* ((1-baseData.weight(bin(i))).*...
            GauNarr + baseData.weight(bin(i)).* GauBro);
    else if i>1
            latProf.flag(i) = 0;
            latProf.Mask(:,:,i) = latProf.Mask(:,:,i-1);
        else
            latProf.Mask(:,:,i) = zeros([length(x), length(x)]);
        end
    end
end
disp('')