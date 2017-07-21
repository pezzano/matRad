vec = [];
for i = 1:stf4.numOfRays
    vec = [vec stf4.ray(i).energy];
end

univec = unique(vec);

for i = 1:length(univec)
    strct.rayPerEnergy{i,1} = univec(i);
    strct.rayPerEnergy{i,2} = [];
    for j = 1:stf4.numOfRays
        for k =  1:length(stf4.ray(j).energy)
            if any(stf4.ray(j).energy(k)==univec(i)) 
                strct.rayPerEnergy{i,2} = [strct.rayPerEnergy{i,2}; stf4.ray(j).rayPos_bev];
            end
        end
    end
end

for i=1:stf4.numOfRays
    strct.rayPos_bev(i,:) = stf4.ray(i).rayPos_bev;
end

stf_new = stf4;
stf_new.rayPos_bev = strct.rayPos_bev;
stf_new.rayPerEnergy = strct.rayPerEnergy;

