ylimRange = [-30,120];
xlimRange = [-.3, .4];
lineWidth = 0.75;
rstar_dot_size = 8; %it was 4 -PP
fastCol = [136 87 164]./255;
regCol = [0 166 81]./255;
alphashade = 0.2;
smoothingKernel = 10;
%%
id = 39;


 ni = goodAndSafe_units_idx(id);          %% ni is index in the sorted unit arrays
       
       unitRasterData = load(fullfile(respRasterDir,...
            sprintf('%s_NPX_unit%d_resps_raster_data.mat',...
            mouseName, ni)));
       
       neuronInhRaster = unitRasterData.raster_data;
       
       neuronRtypeRaster = neuronInhRaster(...
           unitRasterData.raster_labels.labels == 2 &...
           unitRasterData.raster_labels.xpeakCondLabel == 0, :);      
       neuronStypeRaster = neuronInhRaster(...
           unitRasterData.raster_labels.labels == 1 &...
           unitRasterData.raster_labels.xpeakCondLabel == 0, :);      
       neuron_firstsniff_Raster                         = neuronInhRaster(flagX_included_0, :);
       neuron_secondormoresniff_Raster                  = neuronInhRaster(secondormoresniffX_included_0, :);
       neuron_lastregular_Raster                        = neuronInhRaster(lastregular_included_0, :);
       neuron_secondlastregular_Raster                  = neuronInhRaster(secondlastregular_included_0, :);
       neuron_allRbetweenRX_Raster                      = neuronInhRaster(allRbetweenRX_included_0, :);
       
       
       
       
       
       
       %%
figure

slowRast2plot = tableAdd.secondlastregular_Raster{id};
slowPSTH2plot = tableAdd.allRbetweenR{id};
slowTag = 'rRr_inPSTH';


firstsniff = gen_fx_gsmooth((neuron_firstsniff_Raster*1000-...
                    unitsBaselineFR(ni)), smoothingKernel);     
rRr = gen_fx_gsmooth((neuron_allRbetweenRX_Raster*1000-...
                    unitsBaselineFR(ni)), smoothingKernel);            
       


    
    
    subplot(2,1,1)
    hold on
    rasterHeight = 0;
    [row_slow,col_slow] =...
        find(tableAdd.firstsniff_Raster{id} == 1);
    plot(col_slow, row_slow+rasterHeight,'.', 'Color',...
        fastCol, 'MarkerSize', rstar_dot_size);
    rasterHeight = rasterHeight+...
        size(tableAdd.firstsniff_Raster{id}, 1);
    [row_slow,col_slow] =...
        find(slowRast2plot == 1);
    plot(col_slow, row_slow+rasterHeight,'.', 'Color',...
        regCol, 'MarkerSize', rstar_dot_size);
    rasterHeight = rasterHeight+...
        size(slowRast2plot, 1);
    xlim(xlimRange*1000+601)
    ylim([-1,rasterHeight-1])
    yticks([])
    box off
    title(sprintf('%d - %s', id, slowTag))
    
    
    subplot(2,1,2)
    hold on
%     plot(PITH_time, tableAdd.firstsniff{id} +tableAdd.unitsBaselineFR(id),...
%         'Color',fastCol,...
%         'linewidth',lineWidth);
%     plot(PITH_time, slowPSTH2plot + tableAdd.unitsBaselineFR(id),...
%         'Color',regCol,...
%         'linewidth',lineWidth);
    stdshade_modified(firstsniff, alphashade, fastCol ,PITH_time,[],[], 2)
    stdshade_modified(rRr, alphashade, regCol ,PITH_time,[],[], 2)
    
    
    xlim(xlimRange)
    xticks([-.2,0,.2,.4])
    ylabel('Firing rate (Hz)')
    box off
    