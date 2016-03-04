function [Points, distPC2PC, distPC2PCwrapped] = CreateRandomNetwork(cellNo, maxSeperationDistance)

initPointsNo = cellNo;
DistTolerance = 2;
Points = rand(initPointsNo,3);
% Points=sortrows(Points);
%Calculate distances
mDist = zeros((initPointsNo+1)^2,1);
for i=1:initPointsNo
    for j=1:initPointsNo
            %Calculate all distances
            mDist((i-1)*(initPointsNo)+j) = mDist((i-1)*(initPointsNo)+j) + sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2);
    end
end
% Find max value
scrDist = max(mDist) ;
redoDist = (scrDist < maxSeperationDistance+DistTolerance) || (scrDist > maxSeperationDistance-DistTolerance);

%Reshape to max Sepearition Distance 
if redoDist
    DistFactor = scrDist / maxSeperationDistance;
    Points = Points / DistFactor;
end

cubeMax = max(Points(:));
% cubeMin = min(Points);
% cbnIdx = combntns(1:9,2);
cbnOffset = [-cubeMax, 0, cubeMax];
[X,Y,Z] = meshgrid(1:3,1:3,1:3);
cbnOffsets = cbnOffset([X(:),Y(:),Z(:)]);

start = 2;
for i=1:initPointsNo
    for j=start:initPointsNo
        distPC2PC(i,j) = sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2); 
        tmp = zeros(1,length(cbnOffsets));
        for k=1:length(cbnOffsets)
            tmp(k) = sqrt((Points(i,1)-(Points(j,1)+cbnOffsets(k,1)))^2 +...
                (Points(i,2)-(Points(j,2)+cbnOffsets(k,2)))^2 +...
                (Points(i,3)-(Points(j,3)+cbnOffsets(k,3)))^2);
        end
        distPC2PCwrapped(i,j) = min(tmp); 
    end
    start = start+1;
end

% pad with zeros the last row:
distPC2PC(initPointsNo,:) = 0;
distPC2PCwrapped(initPointsNo,:) = 0;

return;