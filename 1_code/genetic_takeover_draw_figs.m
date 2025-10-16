%% Selection of DNA genomes over RNA genomes in early evolution by a stoichiometric selection advantage

% Written in MATLAB R2025a with:
% - Statistics and Machine Learning Toolbox
% - Bioinformatics Toolbox
% - Optimization Toolbox
% - Global Optimization Toolbox
% - Symbolic Math toolbox
% - Parallel Computing Toolbox

% Rosemary Yu. 2025-10-15.


%% load results
close all
clear
clc

load("2_results\ribocell_size_1.mat")
load("2_results\ribocell_size_2.mat")
load("2_results\ribocell_size_3.mat")
load("2_results\kpol_kribo_ratio.mat")

%% violin plots from res1 and res2:

rid = any(res1(:,9) > 0, 2);
s1.rna = res1(rid,6);
rid = any(res1(:,10) > 0, 2);
s1.hyb = res1(rid,7);
rid = any(res1(:,11) > 0, 2);
s1.dna = res1(rid,8);

rid = any(res2(:,9) > 0, 2);
s2.rna = res2(rid,6);
rid = any(res2(:,10) > 0, 2);
s2.hyb = res2(rid,7);
rid = any(res2(:,11) > 0, 2);
s2.dna = res2(rid,8);

fig1 = figure;
violinplot(s1)
ylim([0 9])
ylabel('min volume (x10^{-5} um^3)')
xticks([1 2 3])
xticklabels({'RNA','hybrid','DNA'})
title('ribo and pol only')

fig2 = figure;
violinplot(s2)
ylim([0 14])
ylabel('min volume (x10^{-5} um^3)')
xticks([1 2 3])
xticklabels({'RNA','hybrid','DNA'})
title('with AARS & tRNA')

%saveas(fig1, '3_figures\ribo_pol_only.png')
%saveas(fig2, '3_figures\with_AARS_tRNA.png')


%% composite violin plots & histogram from res3:
close all

fig3 = figure;
%tiledlayout(height(res3), width(res3)) 
tiledlayout(3, width(res3)) % only plot the middle 3 ratios

fP = ["10k" "25k" "100k" "250k" "1M" "2.5M"];
fR = ["2" "1" "1/2"]; % only plot the middle 3 ratios

for i = 3:5 % fR/fP 
    for j = 1:width(res3) % fP size 
        res = res3{i,j};
        
        rid = any(res(:,9) > 0, 2);
        s.rna = res(rid,6);
        rid = any(res(:,10) > 0, 2);
        s.hyb = res(rid,7);
        rid = any(res(:,11) > 0, 2);
        s.dna = res(rid,8);

        if isempty(s.rna)
            s.rna=0;
        end
        if isempty(s.hyb)
            s.hyb=0;
        end
        if isempty(s.dna)
            s.dna=0;
        end

        nexttile
        violinplot(s)
        ylim([0 Inf])
        xlim([0.5 3.5])
        ylabel('min volume (x10^{-5} um^3)')
        xticks([1 2 3])
        xticklabels({'RNA','hybrid','DNA'})
        title("fP/fR = "+ fR(i-2) + "; fP size = " + fP(j))

    end
end

%saveas(fig3, '3_figures\with_enzyme_ribozyme_violins.png')

%% heatmaps from res3:
close all

h_hyb = zeros(height(res3),width(res3));
h_dna = zeros(height(res3),width(res3));

for i = 1:height(res3)
    for j = 1:width(res3)
        res = res3{i,j};
        rid = any(sum(res(:,9:11),2) > 0,2);
        N = res(rid,9:11);
        
        minN1 = min(N(N(:,1)>0,1));
        minN2 = min(N(N(:,2)>0,2));
        minN3 = min(N(N(:,3)>0,3));

        if isempty(minN1)
            minN1 = 1e5; %arbitrary large number
        end
        if isempty(minN2)
            minN2 = 1e5; %arbitrary large number
        end
        if isempty(minN3)
            minN3 = 1e5; %arbitrary large number
        end

        if minN2 < minN1 && minN2 < minN3
            % hybrid advantage
            N2 = N(N(:,2)>0,2);
            h_hyb(i,j) = height(find(N2 == minN2)) / height(N2);
        elseif minN3 < minN1 && minN3 < minN2
            % DNA advantage
            N3 = N(N(:,3)>0,3);
            h_dna(i,j) = height(find(N3 == minN3)) / height(N3);
        end
    end
end

h_hyb = round(h_hyb,2)*100; %in percentage
h_dna = round(h_dna,2)*100; %in percentage


fig4 = figure;
n = 50;
c = [ones(n,1) linspace(1,0.8,n)' linspace(1,0,n)'];  % colormap from white [1 1 1] to light orange [1 0.8 0]
h = heatmap(h_hyb,'ColorLimits',[0 75], 'Colormap', c);
grid off
colorbar off
xlabel("fP size")
ylabel("fR/fP (\pm10%)")
title("hybrid advantage")
h.XDisplayLabels = ["10k" "25k" "100k" "250k" "1M" "2.5M"];
h.YDisplayLabels = ["10" "5" "2" "1" "1/2" "1/5" "1/10"];


fig5 = figure;
n = 50;
c = [linspace(1,0,n)' linspace(1,0.6,n)' linspace(1,0,n)'];  % colormap from white [1 1 1] to green [0 0.6 0]
h = heatmap(h_dna,'ColorLimits',[0 100], 'Colormap', c); 
grid off
colorbar off
xlabel("fP size")
ylabel("fR/fP (\pm10%)")
title("DNA advantage")
h.XDisplayLabels = ["10k" "25k" "100k" "250k" "1M" "2.5M"];
h.YDisplayLabels = ["10" "5" "2" "1" "1/2" "1/5" "1/10"];

%saveas(fig4, '3_figures\with_enzyme_ribozyme_hyb.png')
%saveas(fig5, '3_figures\with_enzyme_ribozyme_dna.png')


%% heatmaps for res4_all:

close all

kpkr = ["100" "20" "3" "1" "1/3" "1/20" "1/100"];

fig6 = figure;
tiledlayout(2,7, 'TileIndexing','columnmajor');
n = 50;
c1 = [ones(n,1) linspace(1,0.8,n)' linspace(1,0,n)'];  % colormap from white [1 1 1] to light orange [1 0.8 0]
c2 = [linspace(1,0,n)' linspace(1,0.6,n)' linspace(1,0,n)'];  % colormap from white [1 1 1] to green [0 0.6 0]

for p = 1:7
    s = res4_all(:,:,p);
    
    h_hyb = zeros(height(s),width(s));
    h_dna = zeros(height(s),width(s));
    
    for i = 1:height(s)
        for j = 1:width(s)
            res = s{i,j};
            rid = any(sum(res(:,9:11),2) > 0,2);
            N = res(rid,9:11);
            
            minN1 = min(N(N(:,1)>0,1));
            minN2 = min(N(N(:,2)>0,2));
            minN3 = min(N(N(:,3)>0,3));
    
            if isempty(minN1)
                minN1 = 1e5; %arbitrary large number
            end
            if isempty(minN2)
                minN2 = 1e5; %arbitrary large number
            end
            if isempty(minN3)
                minN3 = 1e5; %arbitrary large number
            end
    
            if minN2 < minN1 && minN2 < minN3
                % hybrid advantage
                N2 = N(N(:,2)>0,2);
                h_hyb(i,j) = height(find(N2 == minN2)) / height(N2);
            elseif minN3 < minN1 && minN3 < minN2
                % DNA advantage
                N3 = N(N(:,3)>0,3);
                h_dna(i,j) = height(find(N3 == minN3)) / height(N3);
            end
        end
    end
    
    h_hyb = round(h_hyb,2)*100; %in percentage
    h_dna = round(h_dna,2)*100; %in percentage
    
    nexttile
    h = heatmap(h_hyb,'ColorLimits',[0 75], 'Colormap', c1); %orange colormap
    grid off
    colorbar off
    xlabel("fP size")
    ylabel("fR/fP (\pm10%)")
    title("kpol/kribo = "+ kpkr(p) + ", hybrid adv")
    h.XDisplayLabels = ["250k" "1M" "2.5M"];
    %h.YDisplayLabels = ["2" "1" "1/2"];
    h.YDisplayLabels = ["10" "5" "2" "1" "1/2" "1/5" "1/10"];
    
    nexttile
    h = heatmap(h_dna,'ColorLimits',[0 100], 'Colormap', c2); %green colormap
    grid off
    colorbar off
    xlabel("fP size")
    ylabel("fR/fP (\pm10%)")
    title("kpol/kribo = "+ kpkr(p) + ", DNA adv")
    h.XDisplayLabels = ["250k" "1M" "2.5M"];
    %h.YDisplayLabels = ["2" "1" "1/2"];
    h.YDisplayLabels = ["10" "5" "2" "1" "1/2" "1/5" "1/10"];

end

%saveas(fig6, '3_figures\kpol_kribo_ratio.png')