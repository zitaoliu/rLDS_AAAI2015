clc;
clear;

rvals = [255 13 63 13]./255;
gvals = [0 123 255 235]./255;
bvals = [0 255 120 255]./255;

load rLdsRank_Synthetic_models.mat
rLdsRank_Synthetic_models = models;

load rLdsGroupLasso_Synthetic_models.mat
rLdsGroupLasso_Synthetic_models = models;

config = initSynthetic();

figure
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontname', 'Times New Roman')

hiddenStates = config.hiddenStates;
fontSize = 24;

for i = 1:numel(hiddenStates)
    
    A = rLdsRank_Synthetic_models{i, 1}.A;
    sv_rank = svd(A);
    sv_rank = sv_rank./max(sv_rank);
    
    A = rLdsGroupLasso_Synthetic_models{i, 1}.A;
    sv_gl = svd(A);
    sv_gl = sv_gl./max(sv_gl);
    
    svs = [sv_gl sv_rank];
    
    subplot(1, 3, i)
    hold on;
    for k = 1:2
        plot(svs(:, k), '.-', 'LineWidth', 2, 'MarkerSize', 30, ...
            'Color', [rvals(k), gvals(k),  bvals(k)]);
    end
    hold off;
    
    legend({'rLDS_{G}', 'rLDS_{R}'})
    
    xlabel('Singular value indices')
    ylabel('Normalized singluar values')
    title(['# of states: ', num2str(hiddenStates(i))], 'fontsize', fontSize)
    
    set(gca,'FontSize',fontSize)
    
    ylhand = get(gca, 'ylabel'); set(ylhand, 'fontsize', fontSize)
    xlhand = get(gca, 'xlabel'); set(xlhand, 'fontsize', fontSize)
    
end


set(gcf, 'PaperPosition', [0 0 32 8]);
set(gcf, 'PaperSize', [32 8]);
saveas(gca,['synthetic_sv.pdf']);