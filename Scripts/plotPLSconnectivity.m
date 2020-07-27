function [ ] = plotPLSconnectivity(resultRotated,nLV,parcel_indx,thrsh,inds,network_labels,lowTriagDataIndx,ax_limits)
%% displayPLSconnectivity - 
%   Parameters :
%       resultRotated - structured array of the PLS output. The rows in
%           resultRotated.boot_result.compare_u(:)the compared salience or
%           compared brain scores which this function uses as edge-weights in
%           the correlation plot. 
%       nLV -  number of latent variable 
%       parcel_indx - (1 x Nrois) array,index array of ROIs in correlation
%       thrsh - threshold to apply to PLS connectivity matrix.These 
%            operates like z-scores and can be thresholded similarly.
%       inds - array index used to partition matrix by network assignment
%       network - (N,) network assignment for each ROI 
%       lowTriagDataIndx - index of lower triangle of correlation matrix
%       thrsh - threshold to apply to bootstrap ratio 
%       ax_limits - lower and upper limits for the colorbar axes

% NOTE: the colormap used here (ColorBrewer) is not a built-in Matlab 
% colormap. If you do not have it the code will use the deulat (perula).
%% Display bootstrap connectivity values
dataForPlot = struct([]);
for cc = 1 : nLV
    temp1=resultRotated.boot_result.compare_u(:,cc); % second value represents LV of interest --> (:,1) = LV1, (:,2) = LV2, (:,3) = LV3
    %compare_u represents the bootstrap ratio values
    dataForPlot{cc} = zeros(length(parcel_indx));
    dataForPlot{cc} = zeros(length(parcel_indx));
    dataForPlot{cc}(lowTriagDataIndx)= temp1;
    dataForPlot{cc}(abs(dataForPlot{cc})<thrsh)=0; % threshold
    dataForPlot{cc} = dataForPlot{cc} + tril(dataForPlot{cc},-1).'; %copying lower triangle to upper over diagonal
    imagesc(dataForPlot{cc});
    xL = get(gca, 'XLim'); 
    yL = get(gca, 'YLim');
    for n = 2:length(inds) %this does NOT work for matlab versions < 2019.FIX!!
        line([inds(n) inds(n)],yL, 'LineWidth', 2, 'Color', 'k');
        line(xL,[inds(n) inds(n)], 'LineWidth', 2, 'Color', 'k');
    end 
    set(gca,'Layer','top','YTick',inds,'XTick',inds,'YTickLabel',...
        network_labels,'XTickLabel',network_labels);
    xtickangle(45);
    caxis(ax_limits);
    colormap(brewermap([],'*RdBu'));
    colorbar;
    plot_title='PLS Correlation Matrix ';
    title([plot_title, 'of LV: ', num2str(cc)]);
    %Save plot
    filename = sprintf('LV_matrix_%d.png', cc) ;
    saveas(gca, filename) 
end
end

