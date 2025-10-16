%% Evolutionary selection of DNA genomes over RNA genomes by a stoichiometric selection advantage

% Written in MATLAB R2025a with:
% - Statistics and Machine Learning Toolbox
% - Bioinformatics Toolbox
% - Optimization Toolbox
% - Global Optimization Toolbox
% - Symbolic Math toolbox
% - Parallel Computing Toolbox

% Rosemary Yu. 2025-10-15.



%% ribocell size 1. ribo & pol only
% ribocell size 2. with AARS and tRNA
% run time ca. 5 min with 6 parallel workers

clear
clc
close all

tic

kribo = 12;
kpol = kribo*3;
intNflag = 1;

nrep = 500;
res1 = zeros(nrep,19);
res2 = zeros(nrep,19);

parfor i = 1:nrep
    disp(i)
    % ribo & pol only
    nuc_rRNA = randi([4000 4500]); 
    nuc_RP   = randi([6000 9000]); 
    nuc_pol  = randi([1500 1800]); 
    nucCount = [nuc_rRNA nuc_RP nuc_pol];

    [minV1, Ntot1, Nmol1] = rnaGenome(kribo, kpol, intNflag, nucCount);
    [minV2, Ntot2, Nmol2] = hybridGenome(kribo, kpol, intNflag, nucCount);
    [minV3, Ntot3, Nmol3] = dnaGenome(kribo, kpol, intNflag, nucCount);
    res1(i,1:19) = [nucCount 0 0 minV1 minV2 minV3 Ntot1 Ntot2 Ntot3 Nmol1(1:2)' Nmol2(1:3)' Nmol3(1:3)']; 
   
    % add AARS and tRNA
    nuc_AARS = randi([5000 6000]);
    nuc_tRNA = randi([3000 4000]);
    nucCount = [nuc_rRNA nuc_RP nuc_pol nuc_AARS nuc_tRNA];

    [minV1, Ntot1, Nmol1] = rnaGenome(kribo, kpol, intNflag, nucCount);
    [minV2, Ntot2, Nmol2] = hybridGenome(kribo, kpol, intNflag, nucCount);
    [minV3, Ntot3, Nmol3] = dnaGenome(kribo, kpol, intNflag, nucCount);
    res2(i,1:19) = [nucCount minV1 minV2 minV3 Ntot1 Ntot2 Ntot3 Nmol1(1:2)' Nmol2(1:3)' Nmol3(1:3)'];

end

toc

% write output
currentPath = pwd;
cd('2_results')
save('ribocell_size_1.mat', 'res1' );
save('ribocell_size_2.mat', 'res2' );

col_names = {'nuc_rRNA','nuc_RP','nuc_pol','nuc_AARS','nuc_tRNA', ...
    'minV_R','minV_H','minV_G','Ntotal_R', 'Ntotal_H', 'Ntotal_D',...
    'ribo_R','pol_R','ribo_H', 'pol1_H', 'pol2_H', 'ribo_D', 'pol1_D', 'pol2_D'};
res1_table = array2table(res1,'VariableNames',col_names);
res2_table = array2table(res2,'VariableNames',col_names);
writetable (res1_table, 'Table_S1.csv');
writetable (res2_table, 'Table_S2.csv');
cd(currentPath)


%% ribocell size 3. with functional proteins (enzymes) and RNA (ribozymes)
% runtime ca. 1 h with 6 parallel workers

clear
clc
close all

kribo = 12;
kpol = kribo*3;
intNflag = 1;

% scan range:
fprot_range = [1e4 2.5e4 1e5 2.5e5 1e6 2.5e6]; % 2.5e6 = latest estimate LUCA genome size
fR_over_fP_range = [10 5 2 1 1/2 1/5 1/10];

max_i = length(fprot_range);
max_j = length(fR_over_fP_range);
res3 = cell(length(fR_over_fP_range), length(fprot_range)); %each column a fixed fP, each row a fR/fP ratio
nrep = 50; 

tic

parfor i = 1:max_i
    for j = 1:max_j
        res_ij = zeros(nrep,19);
    
        for k = 1:nrep
            disp([i j k])
    
            nuc_rRNA = randi([4000 4500]); 
            nuc_RP   = randi([6000 9000]); 
            nuc_pol  = randi([1500 1800]); 
            nuc_fprot = fprot_range(i); 
            fR_over_fP = randr([fR_over_fP_range(j)*0.9 fR_over_fP_range(j)*1.1]);
            nuc_fRNA = nuc_fprot * fR_over_fP;
            nucCount = [nuc_rRNA nuc_RP nuc_pol nuc_fprot nuc_fRNA];
    
            [minV1, Ntot1, Nmol1] = rnaGenome(kribo, kpol, intNflag, nucCount);
            [minV2, Ntot2, Nmol2] = hybridGenome(kribo, kpol, intNflag, nucCount);
            [minV3, Ntot3, Nmol3] = dnaGenome(kribo, kpol, intNflag, nucCount);
            res_ij(k,1:19) = [nucCount minV1 minV2 minV3 Ntot1 Ntot2 Ntot3 Nmol1(1:2)' Nmol2(1:3)' Nmol3(1:3)'];
        end
        res3{j,i} = res_ij;
    end
end

toc

% write output
currentPath = pwd;
cd('2_results')
save('ribocell_size_3.mat', 'res3' );
cd(currentPath)

% nb. no csv saved (this would write to 6x7=42 csv's)


%% kpol to kribo ratio 
% runtime ca. 30 min with 6 parallel workers

clear
clc
close all

kribo = 12;
intNflag = 1;

% scan range:
kpol_kribo_range = [100 20 1 1/3 1/20 1/100]; % 3 is already done
fprot_range = [2.5e5 1e6 2.5e6]; 
fR_over_fP_range = [10 5 2 1 1/2 1/5 1/10];

max_p = length(kpol_kribo_range);
max_i = length(fprot_range);
max_j = length(fR_over_fP_range);
res4 = cell(length(fR_over_fP_range), length(fprot_range)); 
%every row is an fR/fP ratio, every column is an fP value
%every page is a kpol/kribo ratio 

nrep = 50; 

tic

parfor p = 1:max_p %kpol_kribo_range
    for i = 1:max_i %fprot_range
        for j = 1:max_j %fR_over_fP_range
            res_ij = zeros(nrep,19);
        
            for k = 1:nrep
                disp([p i j k])
        
                nuc_rRNA = randi([4000 4500]); 
                nuc_RP   = randi([6000 9000]); 
                nuc_pol  = randi([1500 1800]); 
                nuc_fprot = fprot_range(i); 
                fR_over_fP = randr([fR_over_fP_range(j)*0.9 fR_over_fP_range(j)*1.1]);
                nuc_fRNA = nuc_fprot * fR_over_fP;
                nucCount = [nuc_rRNA nuc_RP nuc_pol nuc_fprot nuc_fRNA];

                kpol = kribo*kpol_kribo_range(p);
        
                [minV1, Ntot1, Nmol1] = rnaGenome(kribo, kpol, intNflag, nucCount);
                [minV2, Ntot2, Nmol2] = hybridGenome(kribo, kpol, intNflag, nucCount);
                [minV3, Ntot3, Nmol3] = dnaGenome(kribo, kpol, intNflag, nucCount);
                res_ij(k,1:19) = [nucCount minV1 minV2 minV3 Ntot1 Ntot2 Ntot3 Nmol1(1:2)' Nmol2(1:3)' Nmol3(1:3)'];
            end
            res4{j,i,p} = res_ij; 
            %every row is an fR/fP ratio, every column is an fP value
            %every page is a kpol/kribo ratio 
        end
    end
end

toc

% add res3 into it as well (kpol/kribo = 3)
load("2_results\ribocell_size_3.mat")
%res3_trim = res3(3:5,4:6);
res3_trim = res3(:,4:6);

res4_all = res4(:,:,1:2); 
res4_all(:,:,3) = res3_trim; %kp/kr = 3
res4_all(:,:,4:7) = res4(:,:,3:6); 

% write output
currentPath = pwd;
cd('2_results')
save('kpol_kribo_ratio.mat', 'res4_all');
cd(currentPath)

disp('all done')

% nb. no csv saved (this would write to 7x21=147 csv's)



%% local

function y = randr(x)
    % outputs a random non-integer (4 decimal places) from a range given by x = [a b]
    a = x(1);
    b = x(2);
    y = a + (b-a).*rand;

end
