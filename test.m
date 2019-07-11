for i = 1:3
   
    fig = figure;
    x = trials(50+i).Data.EndEffPos_FromJA;
    v = trials(50+i).Data.xdot;
    va = trials(50+i).Data.va;
    
    % find when v crosses 0.02 or -0.02 for the last time before the max v
    [M, I] = max(va);
    inds = abs(va(1:I)) < 0.005;
    firstMoveInd = find(inds, 1, 'last');
    
    inds = abs(va) > 0.001;
    lastMoveInd = find(inds, 1, 'last');
    
    plot(x);
    ax = fig.Children;
    hold on
    plot(v);


    
    ylims = ax.YLim;
    plot(va);    
    % Make line to mark where x has risen 5 percent
    demark_line = zeros(size(x));
%     ind = GetFivePercentRiseInd(x);
    demark_line(1:firstMoveInd) = ylims(1) - 1;
    demark_line(firstMoveInd+1:lastMoveInd) = ylims(2) + 1;
    demark_line(lastMoveInd+1:end) = ylims(1) - 1;
    
    
    plot(demark_line);
    legend('x', 'v', 'va', 'demark')
    ax.YLim = ylims;
    hold off
    
    
    
end