function [Tidal_normal,idx_Tidal_normal] = FxOptimal_Lcurve(curve)
    
%     curve = [8.4 8.3 5.4 5.3 4.8 4.7 4.6 4.51 4.5 4.39 3.8 3.7 3.2 2.8 2.1];
    nPoints = length(curve);
    allCoord = [1:nPoints;curve(:,1)']';
    
%     allCoord = curve(end:-1:1,[1 2]);
    firstPoint = allCoord(1,:);
    lineVec = allCoord(end,:) - firstPoint;
    lineVecN = lineVec / sqrt(sum(lineVec.^2));
    vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
    scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1),2);
    vecFromFirstParallel = scalarProduct * lineVecN;
    vecToLine = vecFromFirst - vecFromFirstParallel;
    distToLine = sqrt(sum(vecToLine.^2,2));
    %     figure; plot(distToLine);
    [~,idxOfBestPoint] = max(distToLine);
    
    if nargin > 1
        figure; plot(curve)
        hold on;
        plot(allCoord(idxOfBestPoint,1), allCoord(idxOfBestPoint,2),'or');
    end
    
    idx_Tidal_normal = allCoord(idxOfBestPoint,1);
    Tidal_normal= allCoord(idxOfBestPoint,2);