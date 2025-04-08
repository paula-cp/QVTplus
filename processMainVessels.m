function locEntry = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT)
    % Initialize parameters
    if length(correspondenceDict.(keyName)) == 1
        locEntry = [correspondenceDict.(keyName), 5];
        return;
    end

    segmentIndices = correspondenceDict.(keyName);
    bestSegment = [];
    highestZ = -inf;
    lowestZ = inf;
    largestY = -inf;
    lowestDist = inf;

    % Loop through each segment index and process it
    for segIdx = segmentIndices'
        segmentPositions = data_struct.branchList(data_struct.branchList(:, 4) == segIdx, 1:3);

        if size(segmentPositions, 1) > 4
            if ismember(keyName, {'RACA', 'LACA'})
                % Find segment with the highest Z value
                maxZ = mean(segmentPositions(:, 3));
                if maxZ > highestZ
                    bestSegment = segIdx;
                    highestZ = maxZ;
                end

            elseif ismember(keyName, {'RPCA', 'LPCA'})
                % % Prioritize based on Y and Z criteria
                % firstPoint = segmentPositions(1, :);
                % lastPoint = segmentPositions(end, :);
                % 
                % if lastPoint(2) + 5 < firstPoint(2) % Check Y condition
                %     meanY = mean(segmentPositions(:, 2));
                %     meanZ = mean(segmentPositions(:, 3));
                %     if meanY + 2 > largestY && meanZ - 15 < lowestZ
                %         bestSegment = segIdx;
                %         largestY = meanY;
                %         lowestZ = meanZ;
                %     end
                % end
                bestSegment = [];
                bestSignalQuality = -inf;
            
                % Merge PCA + PC2 candidates
                segmentIndices = correspondenceDict.(keyName);
                if strcmp(keyName, 'RPCA') && isfield(correspondenceDict, 'RPC2')
                    segmentIndices = [segmentIndices, correspondenceDict.('RPC2')];
                elseif strcmp(keyName, 'LPCA') && isfield(correspondenceDict, 'LPC2')
                    segmentIndices = [segmentIndices, correspondenceDict.('LPC2')];
                end
            
                for segIdx = segmentIndices
                    segmentMask = data_struct.branchList(:, 4) == segIdx;
                    if nnz(segmentMask) < 5
                        continue;
                    end
            
                    flowVals = data_struct.flowPulsatile_val(segmentMask, :);  % [pixels x time]
                    if size(flowVals, 1) > 1
                        corrMatrix = corr(flowVals');  % transpose to [time x pixels], get [pixels x pixels] correlation
                        upperTri = triu(corrMatrix, 1);  % exclude diagonal and lower triangle
                        corrVals = upperTri(upperTri ~= 0);
                        signalStrength = mean(corrVals);  % or median(corrVals)
                    else
                        signalStrength = 0;  % Not enough data to compute correlation
                    end
            
                    if signalStrength > bestSignalQuality
                        bestSegment = segIdx;
                        bestSignalQuality = signalStrength;
                    end
                end
            end
        end
    end

    % Check and update for RPC2/LPC2
    if strcmp(keyName, 'RPCA') && ~isempty(bestSegment) && isfield(correspondenceDict, 'RPC2')
        bestSegment = checkSecondarySegments(bestSegment, correspondenceDict.('RPC2'), data_struct);
    elseif strcmp(keyName, 'LPCA') && ~isempty(bestSegment) && isfield(correspondenceDict, 'LPC2')
        bestSegment = checkSecondarySegments(bestSegment, correspondenceDict.('LPC2'), data_struct);
    end

    % Assign final segment or default to NaN
    if ~isempty(bestSegment)
        locEntry = [bestSegment, 5];
    else
        locEntry = [NaN, NaN];
    end
end

function bestSegment = checkSecondarySegments(bestSegment, secondaryIndices, data_struct)
    largestY2 = -inf;

    for segIdx = secondaryIndices'
        segmentPositions = data_struct.branchList(data_struct.branchList(:, 4) == segIdx, 1:3);

        firstPoint = segmentPositions(1, :);
        lastPoint = segmentPositions(end, :);

        if lastPoint(2) + 5 < firstPoint(2) && (lastPoint(3) + 3) > firstPoint(3) % Conditions on Y and Z
            meanY = mean(segmentPositions(:, 2));
            if meanY > largestY2
                bestSegment = segIdx;
                largestY2 = meanY;
            end
        end
    end
end

