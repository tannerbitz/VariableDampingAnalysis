function res = FilterAndAnalyzeData(trials, TargetDirNum, DampNum, varargin)

    % Parse arguments
    p = inputParser;
    defaultVal = -1;
    addParameter(p, 'SubjectNum', defaultVal);
    
    parse(p, varargin{:});
    
    % Get inputted values to filter by
    subjectNumInput = p.Results.SubjectNum;
    targetDirNumInput = TargetDirNum;
    dampNumInput = DampNum;
    
    % DefaultInds
    groupNumInds = ([trials.GroupNumber] >= 4); % default
    subjectNumInds = logical(ones(1, length(trials)));
    targetDirNumInds = logical(ones(1, length(trials)));
    dampNumInds = logical(ones(1, length(trials)));
    nSamplesInds = ([trials.nSamples] > 2000); % all trials should be at least 2 sec long, thus at least 2000 samples in length
    
    % Refine Inds To Filter With Inputs
    if (subjectNumInput ~= defaultVal)
        subjectNumInds = ([trials.SubjectNumber] == subjectNumInput);
    end
    
    if (targetDirNumInput ~= defaultVal)
        targetDirNumInds = ([trials.TargetDirNum] == targetDirNumInput);
    end
    
    if (dampNumInput ~= defaultVal)
        dampNumInds = ([trials.DampingNumber] == dampNumInput);
    end
    
    % Combine Inds and Filter
    filterInds = (groupNumInds & subjectNumInds & targetDirNumInds & dampNumInds & nSamplesInds);
    trials_set = trials(filterInds);
    
    % Exit function if no trials match these criteria
    if (numel(trials_set) == 0)
        fprintf('No trials with group number >=4, target direction: %d, damping number: %d.\n', ...
                targetDirNum, ...
                dampNum);
        res = [];
        return;
    end
    
    % Find min nSamples in set
    nSamplesSet = [trials_set.nSamples];
    nSamples = min(nSamplesSet);
    
    % Get all end effector positions calc'ed from JA
    x_data = [];
    va_data = [];
    damp_data = [];
    xdot_data = [];
    f_data = [];
    for i = 1:length(trials_set)
        x_data = [x_data, trials_set(i).Data.EndEffPos_FromJA(1:nSamples)];
        va_data = [va_data, trials_set(i).Data.va(1:nSamples)];
        damp_data = [damp_data, trials_set(i).Data.Damping(1:nSamples)]; 
        xdot_data = [xdot_data, trials_set(i).Data.xdot(1:nSamples)];
    end
    
    % Find when movement starts
    firstMoveInds = zeros(size(x_data, 2),1);
    restOfArray = zeros(size(x_data, 2), 1);
    for i = 1:length(trials_set)
        va_temp = trials_set(i).Data.va;
        firstMoveInds(i) = GetFirstMoveInd(va_temp);
        restOfArray(i) = length(va_temp) - firstMoveInds(i);
    end
    
    minFirstMoveInd = min(firstMoveInds);
    minRestOfArray = min(restOfArray);
    x_mm = [];          % movement matching
    va_mm = [];
    damp_mm = [];
    force_mm = [];
    for i = 1:length(trials_set)
        firstInd = firstMoveInds(i);% - minFirstMoveInd + 1;
        lastInd = firstMoveInds(i) + minRestOfArray;
        x_mm = [x_mm, trials_set(i).Data.EndEffPos_FromJA(firstInd:lastInd)];
        va_mm = [va_mm, trials_set(i).Data.va(firstInd:lastInd)];
        damp_mm = [damp_mm, trials_set(i).Data.Damping(firstInd:lastInd)];
        force_mm = [force_mm, trials_set(i).Data.Force(firstInd:lastInd, 1)];
    end

    % Put into results struct
    res = struct;
    res.x_data = x_data;
    res.va_data = va_data;
    res.damp_data = damp_data;
    res.x_mm = x_mm;
    res.va_mm = va_mm;
    res.damp_mm = damp_mm;
    res.force_mm = force_mm;
    
    % Get mean and std
    x_data_mean = mean(x_data, 2);
    x_data_std  = std(x_data, 0, 2);
    va_data_mean    = mean(va_data, 2);
    va_data_std     = std(va_data, 0, 2);
    damp_data_mean  = mean(damp_data, 2);
    damp_data_std   = std(damp_data, 0, 2);
    
    % Put into results struct
    res.x_data_mean = x_data_mean;
    res.x_data_std = x_data_std;
    res.va_data_mean = va_data_mean;
    res.va_data_std = va_data_std;
    res.damp_data_mean = damp_data_mean;
    res.damp_data_std = damp_data_std;
    
    % Filter data - Reject Trials with Violations of +/- 3 std
    % Get upper and lower bounds
    x_data_ub = x_data_mean + 3*x_data_std;
    x_data_lb = x_data_mean - 3*x_data_std;
    va_data_ub = va_data_mean + 3*va_data_std;
    va_data_lb = va_data_mean - 3*va_data_std;
    damp_data_ub = damp_data_mean + 3*damp_data_std;
    damp_data_lb = damp_data_mean - 3*damp_data_std;
    
    % Put into results struct
    res.x_data_lb = x_data_lb;
    res.x_data_ub = x_data_ub;
    res.va_data_lb = va_data_lb;
    res.va_data_ub = va_data_ub;
    res.damp_data_lb = damp_data_lb;
    res.damp_data_ub = damp_data_ub;
    
    % Check for data outside bounds any of x bounds.  If there are any
    % violations, then that trial is excluded.
    trial_violation = zeros(size(x_data,2), 1);
    for i = 1:size(x_data, 2)
        x_data_col = x_data(:,i);
        
        x_ub_violation      = x_data_col    > x_data_ub;
        x_lb_violation      = x_data_col    < x_data_lb;
        
        trial_violation_count = sum(x_ub_violation) + ...
                                sum(x_lb_violation); 
        if (trial_violation_count > 0)
            trial_violation(i) = true;
        else
            trial_violation(i) = false;
        end
    end
    
    % Sort data array into ones with/without violations
    x_data_wv = [];     % with violations
    x_data_wov = [];    % without violations
    va_data_wv = [];
    va_data_wov = [];
    damp_data_wv = [];
    damp_data_wov = [];
    xdot_data_wv = [];
    xdot_data_wov = [];
    
    for i = 1:length(trial_violation)
        if (trial_violation(i) == 0)    % no violations
            x_data_wov(:,end+1) = x_data(:,i);
            va_data_wov(:,end+1) = va_data(:,i);
            damp_data_wov(:,end+1) = damp_data(:,i);
            xdot_data_wov(:,end+1) = xdot_data(:,i);
        else                            % violations
            x_data_wv(:,end+1) = x_data(:,i);
            va_data_wv(:,end+1) = va_data(:,i);
            damp_data_wv(:,end+1) = damp_data(:,i);
            xdot_data_wv(:,end+1) = xdot_data(:,i);
        end
    end
    
    % Put into results struct
    res.x_data_wov = x_data_wov;
    res.x_data_wv = x_data_wv;
    res.va_data_wov = va_data_wov;
    res.va_data_wv = va_data_wv;
    res.damp_data_wov = damp_data_wov;
    res.damp_data_wv = damp_data_wv;
    
    % Recalculate mean and std of data without violations
    x_data_wov_mean = mean(x_data_wov, 2);
    x_data_wov_std  = std(x_data_wov, 0, 2);
    va_data_wov_mean    = mean(va_data_wov, 2);
    va_data_wov_std     = std(va_data_wov, 0, 2);
    damp_data_wov_mean  = mean(damp_data_wov, 2);
    damp_data_wov_std   = std(damp_data_wov, 0, 2);
    
    % Put into res struct
    res.x_data_wov_mean = x_data_wov_mean;
    res.x_data_wov_std = x_data_wov_std;
    res.va_data_wov_mean = va_data_wov_mean;
    res.va_data_wov_std = va_data_wov_std;
    res.damp_data_wov_mean = damp_data_wov_mean;
    res.data_data_wov_std = damp_data_wov_std;
    
    % Recalculate upper and lower bounds from data without violations
    x_data_wov_ub = x_data_wov_mean + 3*x_data_wov_std;
    x_data_wov_lb = x_data_wov_mean - 3*x_data_wov_std;
    va_data_wov_ub = va_data_wov_mean + 3*va_data_wov_std;
    va_data_wov_lb = va_data_wov_mean - 3*va_data_wov_std;
    damp_data_wov_ub = damp_data_wov_mean + 3*damp_data_wov_std;
    damp_data_wov_lb = damp_data_wov_mean - 3*damp_data_wov_std;

    % Put into results struct
    res.x_data_wov_ub = x_data_wov_ub;
    res.x_data_wov_lb = x_data_wov_lb;
    res.va_data_wov_ub = va_data_wov_ub;
    res.va_data_wov_lb = va_data_wov_lb;
    res.damp_data_wov_ub = damp_data_wov_ub;
    res.damp_data_wov_lb = damp_data_wov_lb;
    
    
    % Filter movement matching data
    x_mm_mean = mean(x_mm, 2);
    x_mm_std = std(x_mm, [], 2);
    x_mm_ub = x_mm_mean + 3*x_mm_std;
    x_mm_lb = x_mm_mean - 3*x_mm_std;
    
    x_mm_wov = [];
    x_mm_wv = [];
    va_mm_wov = [];
    va_mm_wv = [];
    damp_mm_wov = [];
    damp_mm_wv = [];
    force_mm_wov = [];
    force_mm_wv = [];
    
    for i = 1:size(x_mm, 2)
        x_temp = x_mm(:,i);
        lb_violations = sum(x_temp < x_mm_lb);
        ub_violations = sum(x_temp > x_mm_ub);
        
        if ( (lb_violations > 0) || (ub_violations > 0 ) )
            x_mm_wv = [x_mm_wv, x_temp];
            va_mm_wv = [va_mm_wv, va_mm(:,i)];
            damp_mm_wv = [damp_mm_wv, damp_mm(:,i)];
            force_mm_wv = [force_mm_wv, force_mm(:,i)];
        else
            x_mm_wov = [x_mm_wov, x_temp];
            va_mm_wov = [va_mm_wov, va_mm(:,i)];
            damp_mm_wov = [damp_mm_wov, damp_mm(:,i)];
            force_mm_wov = [force_mm_wov, force_mm(:,i)];
        end
    end
    
    % Get means, std, bounds of movement matched data without violations 
    x_mm_wov_mean = mean(x_mm_wov, 2);
    x_mm_wov_std = std(x_mm_wov, [], 2);
    x_mm_wov_ub = x_mm_wov_mean + 3*x_mm_wov_std;
    x_mm_wov_lb = x_mm_wov_mean - 3*x_mm_wov_std;
    
    va_mm_wov_mean = mean(va_mm_wov, 2);
    va_mm_wov_std = std(va_mm_wov, [], 2);
    va_mm_wov_ub = va_mm_wov_mean + 3*va_mm_wov_std;
    va_mm_wov_lb = va_mm_wov_mean - 3*va_mm_wov_std;
    
    damp_mm_wov_mean = mean(damp_mm_wov, 2);
    damp_mm_wov_std = std(damp_mm_wov, [], 2);
    damp_mm_wov_ub = damp_mm_wov_mean + 3*damp_mm_wov_std;
    damp_mm_wov_lb = damp_mm_wov_mean - 3*damp_mm_wov_std;

    % Put into struct
    res.x_mm_mean = x_mm_mean;
    res.x_mm_std = x_mm_std;
    res.x_mm_ub = x_mm_ub;
    res.x_mm_lb = x_mm_lb;
    
    res.x_mm_wov = x_mm_wov;
    res.x_mm_wv = x_mm_wv;
    res.va_mm_wov = va_mm_wov;
    res.va_mm_wv = va_mm_wv;
    res.damp_mm_wov = damp_mm_wov;
    res.damp_mm_wv = damp_mm_wv;
    res.force_mm_wov = force_mm_wov;
    res.force_mm_wv = force_mm_wv;
    
    res.x_mm_wov_mean = x_mm_wov_mean;
    res.x_mm_wov_std = x_mm_wov_std;
    res.x_mm_wov_ub = x_mm_wov_ub;
    res.x_mm_wov_lb = x_mm_wov_lb;
    
    res.va_mm_wov_mean = va_mm_wov_mean;
    res.va_mm_wov_std = va_mm_wov_std;
    res.va_mm_wov_ub = va_mm_wov_ub;
    res.va_mm_wov_lb = va_mm_wov_lb;
    
    res.damp_mm_wov_mean = damp_mm_wov_mean;
    res.damp_mm_wov_std = damp_mm_wov_std;
    res.damp_mm_wov_ub = damp_mm_wov_ub;
    res.damp_mm_wov_lb = damp_mm_wov_lb;
    

    % Record where the data that was analyzed came from
    res.SubjectNumbers = unique([trials_set.SubjectNumber]);
    res.TargetDirNum = unique([trials_set.TargetDirNum]);
    res.DampingNum = unique([trials_set.DampingNumber]);
    res.TargetDirText = trials_set(1).TargetDirText;
    res.DampingText = trials_set(1).DampingText;
    res.nSamples = nSamples;
    
    % Calculate Overshoot
    x_data_wov_mean_zero = (mean(x_data_wov_mean(1:30)));
    x_data_wov_mean_ss = (mean(x_data_wov_mean(end-100:end)));
    x_data_wov_mean_range = abs(x_data_wov_mean_ss - x_data_wov_mean_zero);
    [x_data_wov_mean_max, ind] = max(abs(x_data_wov_mean));
    x_data_wov_mean_max = x_data_wov_mean_max*sign(x_data_wov_mean(ind));
    
    res.PercentOvershoot_Mean = (abs(x_data_wov_mean_max - x_data_wov_mean_zero)/x_data_wov_mean_range - 1)*100;
    
    x_data_wov_zero = mean(x_data_wov(1:30,:), 1);
    x_data_wov_ss = mean(x_data_wov(end-100:end,:), 1);
    x_data_wov_range = abs(x_data_wov_ss - x_data_wov_zero);
    [M, I] = max(abs(x_data_wov));
    x_data_wov_max = zeros(size(M));
    for i = 1:length(x_data_wov_max)
        x_data_wov_max(i) = x_data_wov(I(i), i);
    end
    res.PercentOvershoot_data_wov = (abs(x_data_wov_max - x_data_wov_zero)./x_data_wov_range - 1)*100;
    res.PercentOvershoot_data_wov_mean = mean(res.PercentOvershoot_data_wov);
    res.PercentOvershoot_data_wov_std = std(res.PercentOvershoot_data_wov);
    
    % Calculate max velocity
    [M, I] = max(abs(xdot_data_wov),[], 1);
    xdot_data_wov_max = zeros(size(M));
    for i = 1:length(M)
        xdot_data_wov_max(i) = xdot_data_wov(I(i), i);
    end
    res.XdotMax_data_wov = xdot_data_wov_max;
    res.XdotMax_data_wov_mean = mean(abs(xdot_data_wov_max));
    res.XdotMax_data_wov_std = std(abs(xdot_data_wov_max));
    
end