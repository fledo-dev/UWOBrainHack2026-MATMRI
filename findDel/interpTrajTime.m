function [phs_spha,phs_conc] = interpTrajTime(phs_spha,phs_conc,tdwell_in,delTime,datatime)
    % Interpolates phs_spha and phs_conc along the 1st dimension 
    %    Time delay between input samples is tdwell_in
    %    delTime is a time delay that is effectively added to datatime (or
    %       subtracted from the times of the input signals)
    %    Output is at samples defined by datatime
    %
    %   Time for trajectory input is presumed to start at 0.
    %
    % (c) Corey Baron 2021

    if isempty(phs_spha) 
        if isempty(phs_conc)
            error('at least one of phs_spha or phs_conc must be non-empty!')
        end
        lenTraj = size(phs_conc,1);
    else
        lenTraj = size(phs_spha,1);
    end

    trajTim = -delTime + tdwell_in*(0:lenTraj-1); 
    if ~isempty(phs_spha)
        phs_spha = interp1(trajTim, phs_spha, datatime, 'makima',nan); 
    end
    if ~isempty(phs_conc)
        phs_conc = interp1(trajTim, phs_conc, datatime, 'makima',nan); 
    end
end