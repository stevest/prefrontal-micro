function [Points, distMat, distPV2PCwrapped] = CreateCubeNetworkPV(cubeSize, cellNo,PCsomata )
% Initializes neuronal PV somata in a cube of given dimentions
% arg_1 the value of cube side in ?m.
% arg_2 is the number of cells

initPointsNo = cellNo;
DistTolerance = 2;
maxSeperationDistance=cubeSize;

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

% somata = somata .* cubeSize;

distMat = zeros(size(Points,1),size(PCsomata,1));

for i=1:size(Points,1)
   distMat(i,1:size(PCsomata,1)) = reshape(sqrt(sum((PCsomata(:,1:3)-repmat(Points(i,1:3),[size(PCsomata,1),1])).^2,2) ),1,[]);
end

% export a dimension-wrapped version of distance mat:

cubeMax = max(Points(:));
% cubeMin = min(Points);
% cbnIdx = combntns(1:9,2);
cbnOffset = [-cubeMax, 0, cubeMax];
[X,Y,Z] = meshgrid(1:3,1:3,1:3);
cbnOffsets = cbnOffset([X(:),Y(:),Z(:)]);

for i=1:initPointsNo
    for j=1:size(PCsomata,1)
%         distPV2PC(i,j) = sqrt((Points(i,1)-PCsomata(j,1))^2 + (Points(i,2)-PCsomata(j,2))^2 + (Points(i,3)-PCsomata(j,3))^2); 
        tmp = zeros(1,length(cbnOffsets));
        for k=1:length(cbnOffsets)
            tmp(k) = sqrt((Points(i,1)-(PCsomata(j,1)+cbnOffsets(k,1)))^2 +...
                (Points(i,2)-(PCsomata(j,2)+cbnOffsets(k,2)))^2 +...
                (Points(i,3)-(PCsomata(j,3)+cbnOffsets(k,3)))^2);
        end
        distPV2PCwrapped(i,j) = min(tmp); 
    end
end

return;