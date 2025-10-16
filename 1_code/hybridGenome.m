function [hyb_minV, hyb_Ntot, hyb_Nmol] = hybridGenome(kribo, kpol, intNflag, nucCount, roundToNearest)
    % minV              minimum volume in (e-5 um3)
    % Ntot              number of total molecules
    % Nmol              number molecules of each component, in the order of 
    %                   [x(1...4) y(1 2) z]. See definitions in the Methods section of the paper.
    %
    % kribo             AA per ribosome per unit time
    % kpol              nucleotides per polymerase per unit time
    % intNflag          flag for whether or not to perform rounding & return integers
    %                   for N, to avoid fraction of molecules. 1 = yes, 0 = no
    % nucCount          either a 1 x 3 array containing the nucleotide count of 
    %                   [rRNA total_RP pol], or a 1 x 5 array of
    %                   [rRNA total_RP pol y z]. 
    % roundToNearest    rounding factor. Default 0.2.

    % Rosemary Yu. 2025-10-14.


    % set up
    if nargin < 5
        roundToNearest = 0.2; 
    end

    if length(nucCount) == 3 %[rRNA total_RP pol]
        nucCount = [nucCount 0 0]; %[rRNA total_RP pol y z]
    end
    aaCount = [nucCount(2)/3 nucCount(3)/3 nucCount(4)/3]; % [RP polymerase y]
    aaCount = floor(aaCount); % in case nucCount is not divisible by 3

    % main
    % Create variables & set init points
    x = optimvar("x",10,1,"LowerBound",1);
    C = optimvar("C","LowerBound",0);
    Ndnuc = optimvar("Ndnuc","LowerBound",1);
    Nrnuc = optimvar("Nrnuc","LowerBound",1);
    Naa = optimvar("Naa","LowerBound",1);

    initialPoint.x = ones(size(x));
    initialPoint.C = 1;
    initialPoint.Ndnuc = 1;
    initialPoint.Nrnuc = 1;
    initialPoint.Naa = 1;

    if nucCount(4) == 0
        y = optimvar("y",3,1,"LowerBound",0, "UpperBound",0);
        initialPoint.y = zeros(size(y));
    else
        y = optimvar("y",3,1,"LowerBound",1);
        initialPoint.y = ones(size(y));
    end

    if nucCount(5) == 0
        z = optimvar("z",2,1,"LowerBound",0, "UpperBound",0);
        initialPoint.z = zeros(size(z));
    else
        z = optimvar("z",2,1,"LowerBound",1);
        initialPoint.z = ones(size(z));
    end

    % Create first LP
    prob = optimproblem("ObjectiveSense","min");
    prob.Objective = sum(x) + sum(y) + sum(z);

    % Define constraints
    prob.Constraints.c1 = Naa == C*kribo*x(1);
    prob.Constraints.c2 = Ndnuc == C*kpol*x(2);
    prob.Constraints.c3 = Nrnuc == C*kpol*x(3);
    prob.Constraints.c4 = Naa == aaCount(1)*x(1) + aaCount(2)*( x(2)+x(3) ) + aaCount(3)*y(1);
    prob.Constraints.c5 = Ndnuc == nucCount(1)*x(7) + nucCount(2)*x(8) + nucCount(3)*( x(9)+x(10) ) + ...
                                   nucCount(4)*y(3) + nucCount(5)*z(2);
    prob.Constraints.c6 = Nrnuc == nucCount(1)*x(1) + nucCount(2)*x(4) + nucCount(3)*( x(5)+x(6) ) + ...
                                   nucCount(4)*y(2) + nucCount(5)*z(1);

    % handle warning 
    warning('error', 'MATLAB:nearlySingularMatrix');

    % solve first LP
    try
        [~,min_Ntotal,exitflag1] = solve(prob,initialPoint, options = optimset('Display','off'));
    catch
        exitflag1 = 0;
    end

   % set up second LP
    if exitflag1 > 0
        prob.Constraints.c7 = sum(x) + sum(y) + sum(z) <= min_Ntotal*1.01; % 1% flexibility
        prob.Objective = x(7) + x(8) + x(9) + x(10) + y(3) + z(2); 
        
        % solve second LP
        try
            [sol,~,exitflag2] = solve(prob,initialPoint, options = optimset('Display','off'));
        catch
            exitflag2 = 0;
        end
    end    
     
    % handle warning
    warning('on', 'MATLAB:nearlySingularMatrix');

    % prep output
    if exitflag1 > 0 && exitflag2 > 0
        res = [sol.x; sol.y; sol.z];
        minV = 300*(sol.Ndnuc + sol.Nrnuc) + 140*sol.Naa; %atomic volumes for nucleotides and AAs
        rf = 1/roundToNearest; 
        res = round(res*rf)/rf; 

        if intNflag == 1
            denom = gcd(sym(res));
            res = res/denom;
            minV = minV/denom;
        end

        hyb_minV = double(minV) * 1e-12 / 1e-5; %times 1e-12 to convert to um3, divide by 1e-5 to manage output
        hyb_Ntot = double(sum(res));
        if nucCount(4) == 0 && nucCount(5) == 0
            hyb_Nmol = double(res(1:10));
        else
            hyb_Nmol = double(res);
        end
        
    else
        hyb_minV = 0;
        hyb_Ntot = 0;
        hyb_Nmol = zeros(15,1);
    end


end
