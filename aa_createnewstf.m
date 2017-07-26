vec = [];
for i = 1:stf.numOfRays
    vec = [vec stf.ray(i).energy];
end

univec = unique(vec);

for i = 1:length(univec)
    strct.rayPerEnergy{i,1} = univec(i);
    strct.rayPerEnergy{i,2} = [];
    for j = 1:stf.numOfRays
        for k =  1:length(stf.ray(j).energy)
            if any(stf.ray(j).energy(k)==univec(i)) 
                strct.rayPerEnergy{i,2} = [strct.rayPerEnergy{i,2}; stf.ray(j).rayPos_bev];
            end
        end
    end
end

for i=1:stf.numOfRays
    strct.rayPos_bev(i,:) = stf.ray(i).rayPos_bev;
end

stf_new3 = stf;
stf_new3.rayPos_bev = strct.rayPos_bev;
stf_new3.rayPerEnergy = strct.rayPerEnergy;

