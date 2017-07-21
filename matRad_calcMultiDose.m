function [doseMat] = matRad_calcMultiDose(TracMat,baseData,latProf,V,tot)

conversionFactor = 1.6021766208e-02;

TracMatc = TracMat;
cutt = any(TracMat,1);
TracMatc( :, ~cutt) = [];
cutt2 = any(TracMatc,2);
TracMatc(~cutt2,:) = [];

dimTM = size(TracMatc);

[x,~]=meshgrid(1:dimTM(2),1:dimTM(1));
[x1,~]=meshgrid(1:dimTM(2),1:length(baseData.Z));

if dimTM(2) == 1
    interpData = interp1(repmat(baseData.depths + baseData.offset,[1 dimTM(2)]),...
        repmat(conversionFactor.*baseData.Z,[1 dimTM(2)]),TracMatc);
else
    interpData = griddata(repmat(baseData.depths + baseData.offset,[1 dimTM(2)]),x1,...
        repmat(conversionFactor.*baseData.Z,[1 dimTM(2)]),TracMatc,x);
end

TracMatf = zeros(size(TracMat));
TracMatf(cutt2,cutt) = interpData;
% resul = zeros(size(TracMatf));
% tic
% for i =1:5
%     convMat = TracMatf;
%     convMat(TracMatf <= latProf.section(i,1) | TracMatf > latProf.section(i,2)) = 0;
%     resul = resul + convn(convMat,permute(latProf.Mask(:,:,i),[3,2,1]),'same');
% end
% toc

TracMatf(~ismember(TracMatf,TracMatf(V))) = 0;

mat4conv = TracMatf;
cut1 = any(any(TracMatf,1),2);
k1 = find(cut1); cut1(k1(1):k1(end)) = 1; k1 = find(cut1);
mat4conv( :,:, ~cut1) = [];

cut2 = any(any(mat4conv,2),3);
k2 = find(cut2); cut2(k2(1):k2(end)) = 1; k2 = find(cut2);
mat4conv(~cut2,:,:) = [];

cut3 = any(any(mat4conv,1),3);
k3 = find(cut3); cut3(k3(1):k3(end)) = 1; k3 = find(cut3);
mat4conv(:,~cut3,:) = [];

doseMat = zeros(size(TracMat));

maxInSampl = max(mat4conv,[],1);
for i =1:tot
    convMat = mat4conv;
    convMat(mat4conv > repmat(maxInSampl*i/tot,size(mat4conv,1),1,1) |...
        mat4conv <= repmat(maxInSampl*(i-1)/tot,size(mat4conv,1),1,1)) = 0;
    if sum(convMat) ~= 0
        resul = convn(convMat,permute(latProf.Mask(:,:,i),[3,2,1]));
        
        cut1(k1(1)-floor(size(latProf.Mask(:,:,i),1)/2) : k1(end)+floor(size(latProf.Mask(:,:,i),1)/2)) = 1;
        cut3(k3(1)-floor(size(latProf.Mask(:,:,i),1)/2) : k3(end)+floor(size(latProf.Mask(:,:,i),1)/2)) = 1;
        doseMat(cut2,cut3,cut1) = doseMat(cut2,cut3,cut1) + resul;
    end
end

