function [TracMat,posShiftNorep,finalWeights]=matRad_multipleRayTracing(radDepthMat,rayEn,iso,res,Dx,Dz,W)

Dx = Dx./res(1);
Dz = Dz./res(3);

pos = rayEn{2} + repmat(iso,[size(rayEn{2},1) 1]);
pos = pos ./repmat(res,[size(pos,1) 1]);

temp = bsxfun(@plus,pos(:,1),Dx);
posShift(:,1) = round(reshape(temp',[],1));
temp = bsxfun(@plus,pos(:,3),Dz);
posShift(:,3) = round(reshape(temp',[],1));
W = repmat(W',[size(pos,1) 1]);

% Check for repetitions of samples
[posShiftNorep,ia,ic] = unique(posShift,'rows','stable');
% posShiftRep = posShift(~ismember([1:length(ic)],ia),:);

finalWeights = W(ia);
finalWeights(ic(~ismember([1:length(ic)],ia))) = ...
    finalWeights(ic(~ismember([1:length(ic)],ia))) + W(~ismember([1:length(ic)],ia));

dim = size(radDepthMat);

pos2 = repmat(posShiftNorep,[1 1 dim(2)]);

x = 1:dim(2);

xx = repmat(x,[size(posShiftNorep,1) 1]);

pos2(:,2,:) = xx;

idx = sub2ind(dim, reshape(pos2(:,2,:),1,[]), reshape(pos2(:,1,:),1,[]), reshape(pos2(:,3,:),1,[]));

TracMat = zeros(size(radDepthMat));
TracMat(idx) = radDepthMat(idx);