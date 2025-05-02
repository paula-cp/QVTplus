function [correspondenceDict, LOCs] = generateLOCs(data_struct, correspondenceDict, multiQVT)
    % Generate the LOCs structure based on the data_struct and correspondenceDict
    % Handles ICA, BA cutoff logic, and main vessel LOC identification

    LOCs = struct(); % Initialize the LOCs structure

    %% Step 1: Process ICA and BA vessels
    ica_slice = 0;
    LICA = correspondenceDict.LICA;
    RICA = correspondenceDict.RICA;
    BA = correspondenceDict.BASI;

    % Extract branch information for each vessel
    info_LICA = find_LOCs('extractBranchInfo', data_struct, LICA);
    info_RICA = find_LOCs('extractBranchInfo', data_struct, RICA);
    info_BA = find_LOCs('extractBranchInfo', data_struct, BA);

    % Round Z values
    info_LICA(:, 3) = round(info_LICA(:, 3));
    info_RICA(:, 3) = round(info_RICA(:, 3));
    info_BA(:, 3) = round(info_BA(:, 3));

    % Find unique Z values
    unique_z_LICA = unique(info_LICA(:, 3));
    unique_z_RICA = unique(info_RICA(:, 3));
    unique_z_BA = unique(info_BA(:, 3));

    % Find the largest Z slice with one point for each vessel
    max_single_z = find_LOCs('findMaxSingleZ', unique_z_LICA, unique_z_RICA, unique_z_BA, ...
                             info_LICA, info_RICA, info_BA);

    if max_single_z > -inf
        ica_slice = max_single_z;

        % Extract locations corresponding to the slice
        LICA_LOC = find_LOCs('extractLocation', info_LICA, max_single_z, 3);
        RICA_LOC = find_LOCs('extractLocation', info_RICA, max_single_z, 3);
        BA_LOC = find_LOCs('extractLocation', info_BA, max_single_z, 3);

        % Ensure minimum LOC Z offset of 3
        LICA_LOC = find_LOCs('ensureMinZOffset', info_LICA, LICA_LOC, 4, 2); %si siguen apareciendo vessels cortados, subir el 4 a 5 o 6.
        RICA_LOC = find_LOCs('ensureMinZOffset', info_RICA, RICA_LOC, 4, 2);
        BA_LOC = find_LOCs('ensureMinZOffset', info_BA, BA_LOC, 4, 2);

        % Store the LOCs for ICA and BA
        LOCs.LICA = [LICA_LOC(1, 4), LICA_LOC(1, 5)];
        LOCs.RICA = [RICA_LOC(1, 4), RICA_LOC(1, 5)];
        LOCs.BASI = [BA_LOC(1, 4), BA_LOC(1, 5)];
    end

    %% Step 2: Process other vessels
    vesselLabels = fieldnames(correspondenceDict);
    for i = 1:numel(vesselLabels)
        keyName = vesselLabels{i};

        if strcmp(keyName, 'LICA') || strcmp(keyName, 'RICA') || strcmp(keyName, 'BASI')
            % Already processed ICA and BA
            continue;
        elseif ismember(keyName, {'LPCA', 'RPCA', 'LMCA', 'RMCA', 'RACA', 'LACA'})
            % Process main vessels
            LOCs.(keyName) = processMainVessels(keyName, correspondenceDict, data_struct, multiQVT);
        elseif strcmp(keyName, 'COMM')
            % Special case for COMM vessels
            LOCs.(keyName) = unique(correspondenceDict.COMM);
        end
    end

    %% Step 3: Handle PCA and secondary PCA logic
    if isfield(correspondenceDict, 'RPC2') && isfield(LOCs, 'RPCA')
        if ismember(LOCs.('RPCA')(1), correspondenceDict.('RPC2'))
            correspondenceDict.('RPCA') = [correspondenceDict.('RPCA'); correspondenceDict.('RPC2')];
            correspondenceDict = rmfield(correspondenceDict, 'RPC2');
        else
            correspondenceDict = rmfield(correspondenceDict, 'RPC2');
        end
    end
    
    if isfield(correspondenceDict, 'LPC2') && isfield(LOCs, 'LPCA')
        if ismember(LOCs.('LPCA')(1), correspondenceDict.('LPC2'))
            correspondenceDict.('LPCA') = [correspondenceDict.('LPCA'); correspondenceDict.('LPC2')];
            correspondenceDict = rmfield(correspondenceDict, 'LPC2');
        else
            correspondenceDict = rmfield(correspondenceDict, 'LPC2');
        end
    end

    %% Step 4: Process communicating arteries
    commNames = {'ACOM', 'LPCO', 'RPCO'};
    for i = 1:numel(commNames)
        keyName = commNames{i};

        if ~isfield(correspondenceDict, keyName)
            continue;
        end

        segmentIDs = correspondenceDict.(keyName);

        % Count how many points each segment has
        lengths = arrayfun(@(j) sum(data_struct.branchList(:, 4) == j), segmentIDs);

        % Filter to segments longer than 5 points
        validIdx = find(lengths > 5);

        if isempty(validIdx)
            warning('%s: No valid segment with >5 points found. Skipping.', keyName);
            continue;
        end

        % Choose the shortest valid one
        [~, minIdx] = min(lengths(validIdx));
        selectedBranch = segmentIDs(validIdx(minIdx));

        % Get info for the selected branch
        length = sum(data_struct.branchList(:,4) == selectedBranch);

        LOCs.(keyName) = [selectedBranch, round(length / 2)];
    end


end
