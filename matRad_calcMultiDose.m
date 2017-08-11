function [doseMat] = matRad_calcMultiDose(radDepthPoints,baseData,latProf,V,gantryAngle,weights,posIdx,dim)

conversionFactor = 1.6021766208e-02;

% TracMatc = TracMat;
% cutt = any(TracMat,1);
% TracMatc( :, ~cutt) = [];
% cutt2 = any(TracMatc,2);
% TracMatc(~cutt2,:) = [];
%
% dimTM = size(TracMatc);
%
% [x,~]=meshgrid(1:dimTM(2),1:dimTM(1));
% [x1,~]=meshgrid(1:dimTM(2),1:length(baseData.Z));

% if dimTM(2) == 1
%     interpData = interp1(repmat(baseData.depths + baseData.offset,[1 dimTM(2)]),...
%         repmat(conversionFactor.*baseData.Z,[1 dimTM(2)]),TracMatc);
% else
%     interpData = griddata(repmat(baseData.depths + baseData.offset,[1 dimTM(2)]),x1,...
%         repmat(conversionFactor.*baseData.Z,[1 dimTM(2)]),TracMatc,x);
% end

interpData = interp1(baseData.depths + baseData.offset,...
    conversionFactor.*baseData.Z,radDepthPoints,'linear');

% TracMatRes = zeros(size(TracMat));
% TracMatRes(cutt2,cutt) = interpData;
% resul = zeros(size(TracMatf));
% tic
% for i =1:5
%     convMat = TracMatf;
%     convMat(TracMatf <= latProf.section(i,1) | TracMatf > latProf.section(i,2)) = 0;
%     resul = resul + convn(convMat,permute(latProf.Mask(:,:,i),[3,2,1]),'same');
% end
% toc

% TracMatf(~ismember(TracMatf,TracMatf(V))) = 0;

% TracMatRes = zeros(dim);
% TracMatRes(posIdx) = interpData .* weights;
%
% TracMatf = zeros(dim);
% TracMatf(V) = TracMatRes(V);
% TracMatf(isnan(TracMatf)) = 0;


% for i = 1:length(weights)
%     TracMatf(:,pos(i,1),pos(i,3)) = TracMatf(:,pos(i,1),pos(i,3)).* weights(i);
% end

% mat4conv = zeros(dim);
% mat4conv(posIdx) = interpData;
% Vmat = zeros(dim); Vmat(V) = 1;
% mat4conv = mat4conv .* Vmat;
% cut3 = any(any(mat4conv,1),2);
% k3 = find(cut3); cut3(k3(1):k3(end)) = 1; k3 = find(cut3);
% mat4conv( :,:, ~cut3) = [];
% 
% cut1 = any(any(mat4conv,2),3);
% k1 = find(cut1); cut1(k1(1):k1(end)) = 1; k1 = find(cut1);
% mat4conv(~cut1,:,:) = [];
% 
% cut2 = any(any(mat4conv,1),3);
% k2 = find(cut2); cut2(k2(1):k2(end)) = 1; k2 = find(cut2);
% mat4conv(:,~cut2,:) = [];

doseMat = zeros(dim);

% find a method the rounds all the numbers for its correct value
[urd,~,ic] = unique(round(radDepthPoints,0),'rows','stable');

% changing index in mat4conv system
% [coord(:,1),coord(:,2),coord(:,3)] = ind2sub(dim,posIdx);
% 
% cutCoord1 = coord(:,1)>max(k1) | coord(:,2)>max(k2) | coord(:,3)>max(k3) ;
% coord(cutCoord1 ,:) = [];
% 
% coord(:,1) = coord(:,1) - min(k1)+1;
% coord(:,2) = coord(:,2) - min(k2)+1;
% coord(:,3) = coord(:,3) - min(k3)+1;
% 
% cutCoord2 = coord(:,1)<=0 | coord(:,2)<=0 | coord(:,3)<=0 ;
% coord(cutCoord2 ,:) = [];
% 
% % coord = round(coord);
% posIdxConv = sub2ind(size(mat4conv), coord(:,1),coord(:,2),coord(:,3));
% 
% interpData (cutCoord1) = []; interpData (cutCoord2) = [];
% ic (cutCoord1) = []; ic (cutCoord2) = [];

% maxInSampl = max(mat4conv,[],1);

% nanidx = isnan(interpData);
% radDepthPoints(nanidx) = [];
% extra = ic(nanidx);

% this does the first computation with all zeros...
for i = 1 :size(latProf.Mask,3)
    if isnan(sum(sum(latProf.Mask(:,:,i))))
        continue;
    end
    %     convMat = mat4conv;
    %     convMat(mat4conv > repmat(maxInSampl*i/tot,size(mat4conv,1),1,1) |...
    %         mat4conv <= repmat(maxInSampl*(i-1)/tot,size(mat4conv,1),1,1)) = 0;
    %     convMat = zeros(size(mat4conv));
    %     convMat(posIdxConv(ic == i)) = interpData(ic == i) .* weights(ic == i);
    
        
        mat4conv = zeros(dim);
        mat4conv(posIdx(ic == i)) = interpData(ic == i) .* weights(ic == i);
%         Vmat = zeros(dim); Vmat(V) = 1;
%         mat4conv = mat4conv .* Vmat;
    
    if sum(sum(sum(mat4conv))) ~= 0 
        cut3 = any(any(mat4conv,1),2);
        k3 = find(cut3); cut3(k3(1):k3(end)) = 1; k3 = find(cut3);
        mat4conv( :,:, ~cut3) = [];
        
        cut1 = any(any(mat4conv,2),3);
        k1 = find(cut1); cut1(k1(1):k1(end)) = 1; k1 = find(cut1);
        mat4conv(~cut1,:,:) = [];
        
        cut2 = any(any(mat4conv,1),3);
        k2 = find(cut2); cut2(k2(1):k2(end)) = 1; k2 = find(cut2);
        mat4conv(:,~cut2,:) = [];
        
        
        mask = imrotate(permute(latProf.Mask(:,:,i),[3,2,1]),-gantryAngle,'bicubic');
        
        resul = convn(mat4conv,mask);
        
        cut3(k3(1)-floor(size(mask,3)/2) : k3(end)+floor(size(mask,3)/2)) = 1;
        cut2(k2(1)-floor(size(mask,2)/2) : k2(end)+floor(size(mask,2)/2)) = 1;
        cut1(k1(1)-floor(size(mask,1)/2) : k1(end)+floor(size(mask,1)/2)) = 1;
        doseMat(cut1,cut2,cut3) = doseMat(cut1,cut2,cut3) + resul;
        if i == 80
            disp('')
        end
    end
end

disp('')
