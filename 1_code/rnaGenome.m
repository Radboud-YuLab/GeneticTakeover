function [rna_minV, rna_Ntot, rna_Nmol] = rnaGenome(kribo, kpol, intNflag, nucCount, roundToNearest)
    % minV              minimum volume in (e-5 um3)
    % Ntotal            number of total molecules
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
    x = optimvar("x",7,1,"LowerBound",1);
    C = optimvar("C","LowerBound",0);
    Nnuc = optimvar("Nnuc","LowerBound",1);
    Naa = optimvar("Naa","LowerBound",1);

    initialPoint.x = ones(size(x));
    initialPoint.C = 1;
    initialPoint.Nnuc = 1;
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
    %prob.Constraints.c1 = Naa == C*kribo*x(1)*f; %see next line for f
    prob.Constraints.c1 = Naa == C*kribo*x(1) * ( 0.5*Nnuc - nucCount(1)*x(1) - nucCount(5)*z(1) ) / ( Nnuc - nucCount(1)*x(1) - nucCount(5)*z(1));
    prob.Constraints.c2 = Nnuc == C*kpol*x(2);
    prob.Constraints.c3 = Naa == aaCount(1)*x(1) + aaCount(2)*x(2) + aaCount(3)*y(1);
    prob.Constraints.c4 = Nnuc == nucCount(1)*(x(1)+x(5)) + nucCount(2)*(x(3)+x(6)) + nucCount(3)*(x(4)+x(7)) + ...
                                      nucCount(4)*(y(2)+y(3)) + nucCount(5)*(z(1)+z(2));
    
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
        prob.Constraints.c5 = sum(x) + sum(y) + sum(z) <= min_Ntotal*1.01; % 1% flexibility
        prob.Objective = x(5) + x(6) + x(7) + y(3) + z(2); 
        
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
        minV = 300*sol.Nnuc + 140*sol.Naa; %300 and 140 are atomic volumes for nucleotides and AAs
        rf = 1/roundToNearest; 
        res = round(res*rf)/rf; 

        if intNflag == 1
            denom = gcd(sym(res));
            res = res/denom;
            minV = minV/denom;
        end

        rna_minV = double(minV) * 1e-12 / 1e-5; %times 1e-12 to convert to um3, divide by 1e-5 to manage output
        rna_Ntot = double(sum(res));
        
        if nucCount(4) == 0 && nucCount(5) == 0
            rna_Nmol = double(res(1:7));
        else
            rna_Nmol = double(res);
        end
        
    else
        rna_minV = 0;
        rna_Ntot = 0;
        rna_Nmol = zeros(12,1);
    end

end

