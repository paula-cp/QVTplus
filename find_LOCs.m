function varargout = find_LOCs(funcName, varargin)
    % Dispatcher function
    switch funcName
        case 'extractBranchInfo'
            varargout{1} = extractBranchInfo(varargin{:});
        case 'findMaxSingleZ'
            varargout{1} = findMaxSingleZ(varargin{:});
        case 'extractLocation'
            varargout{1} = extractLocation(varargin{:});
        case 'ensureMinZOffset'
            varargout{1} = ensureMinZOffset(varargin{:});
        otherwise
            error('Unknown function name: %s', funcName);
    end
end

function branchInfo = extractBranchInfo(data_struct, branchLabels)
    branchInfo = [];
    for i = 1:length(branchLabels)
        branchInfo = [branchInfo; data_struct.branchList(data_struct.branchList(:, 4) == branchLabels(i), :)];
    end
end

function maxZ = findMaxSingleZ(z_LICA, z_RICA, z_BA, info_LICA, info_RICA, info_BA)
    maxZ = -inf;
    for z = max([z_LICA; z_RICA; z_BA]):-1:min([z_LICA; z_RICA; z_BA])
        count_LICA = sum(info_LICA(:, 3) == z);
        count_RICA = sum(info_RICA(:, 3) == z);
        count_BA = sum(info_BA(:, 3) == z);

        if count_LICA == 1 && count_RICA == 1 && count_BA == 1
            maxZ = z;
            break;
        end
    end
end

function loc = extractLocation(info, zValue, zColumn)
    loc = info(find(info(:, zColumn) == zValue, 1), :);
end

function loc = ensureMinZOffset(info, loc, minOffset)
    % Ensure minimum Z offset is at least minOffset
    if loc(1, 5) < minOffset
        loc(1, 5) = minOffset;
    end

    % Calculate max Z value for the same branch
    maxZ = max(info(info(:, 4) == loc(4), 5));

    % Adjust if the difference is less than minOffset
    if (maxZ - loc(1, 5)) < minOffset
        loc(1, 5) = maxZ - minOffset;
    end
end
