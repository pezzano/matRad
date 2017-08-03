function [idx,idx_shift] = mR_shift(initIx,dim,source,target,iso, Dx, Dz)
% This function project a point on a certain ray and return the index of
% the projected point in the reference system of the ct cube
%
% The inputs are:
%   - initIx : initial indeces of the points
%   - cubeDim : dimension of the ct cube (i.e. ct.cubeDim)
%   - source : source point of the ray
%   - target : target point of the ray
%   - iso : isocenter coordinates
%   - Dx : displacement on x axis
%   - Dz : displacement on z axis
% 
% The outputs are:
%   - idx : projected indeces
%   - idx_shift : indeces of the ray translated on the secondary sub-beam

% iso = [stf.isoCenter(2),stf.isoCenter(1),stf.isoCenter(3)];
% source = stf(1).sourcePoint_bev + iso;
% target = stf(1).ray(2).targetPoint_bev + iso;
% dim = ct.cubeDim;
% initIx = V(ix);

% shift for isocenter and new ray
target = target + iso;
target(1) = target(1) + Dx;
target(3) = target(3) + Dz;
source = source + iso;
source(1) = source(1) + Dx;
source(3) = source(3) + Dz;

% convert index in coordinates
[coord(:,1),coord(:,2),coord(:,3)] = ind2sub(dim,initIx);

Avec = repmat(source,length(initIx),1);
Bvec = repmat(target,length(initIx),1);

% calculates the distance between the point and nearest one on the ray; and
% the distance between the nearest point on the ray and one extreme of the
% ray
perpDist =  sqrt(sum(cross(coord-Avec,coord-Bvec).^2,2))./sqrt(sum((Bvec-Avec).^2,2));
otherDist = sqrt( sqrt(sum((coord-Bvec).^2,2)).^2 - perpDist.^2 );

% angles in spherical coordinates
theta = asin(sqrt(sum((Avec(:,3)-Bvec(:,3)).^2,2))./sqrt(sum((Avec-Bvec).^2,2)));
phi = acos(sqrt(sum((Avec(:,1)-Bvec(:,1)).^2,2))./sqrt(sum((Avec(:,1:2)-Bvec(:,1:2)).^2,2)));

% add translation to the extreme of the ray, according to spherical coord
D = round( Bvec - [ -otherDist.*cos(theta).*cos(phi) otherDist.*cos(theta).*sin(phi) -otherDist.*sin(theta)]);

% index the found coordinates
idx = round(sub2ind(dim,D(:,1),D(:,2),round(D(:,3)./10)));

% Shift all indeces 
idx_shift = sub2ind(dim, coord(:,1)+Dx, coord(:,2)-( D(:,2)-coord(:,2) ), coord(:,3)+Dz./10 );

