clear
clc
workdir='G:\brainwander_gitea\jobs\2021_caichao_clustering\Clustering_results_final';
PANSS_hier=[workdir,'\PANSS_clustering\20210923_hierarchical_subtype_combat_5fold_CalinskiHarabasz.mat'];
PANSS_kmeans=[workdir,'\PANSS_clustering\20210923_kmeans_subtype_combat_5fold_CalinskiHarabasz.mat'];
GMV_hier=[workdir,'\GMV_clustering\hierarchical_subtype_results_GMV_5fold_CalinskiHarabasz.mat'];
GMV_kmeans=[workdir,'\GMV_clustering\kmeans_subtype_results_GMV_5fold_CalinskiHarabasz.mat'];
datapath={PANSS_hier,PANSS_kmeans,GMV_hier,GMV_kmeans};
dataname={'PANSS_hier','PANSS_kmeans','GMV_hier','GMV_kmeans'};
for d=1:length(datapath)
    load(datapath{d});
    %plot ARI-CH plot
    ARI_CH_labels={};
    if d==1||d==3
        for h1=1:length(hparams.hdist)
            for h2=1:length(hparams.hmethod)
                ARI_CH_labels{h1,h2}=[hparams.hdist{h1},'-',hparams.hmethod{h2}];
            end
        end
        m=size(ARI_CH_hier);
        y=reshape(ARI_CH_hier,m(1),m(2)*m(3));
    else
        for h1=1:length(kparams.kdist)
            for h2=1:length(kparams.kstart)
                ARI_CH_labels{h1,h2}=[kparams.kdist{h1},'-',kparams.kstart{h2}];
            end
        end
        m=size(ARI_CH_kmeans);
        y=reshape(ARI_CH_kmeans,m(1),m(2)*m(3));
    end
    ARI_CH_labels=reshape(ARI_CH_labels,1,m(2)*m(3));
    maxY=max(y);
    [~,id]=max(maxY');
    x=(2:m(1)+1);
    figure(1);
    hold on
    
    for k=1:m(2)*m(3)
        eval(['hCH_',num2str(k),'=plot(x,y(:,k),''LineWidth'',2);']);
    end
    eval(['legend(hCH_',num2str(id),',ARI_CH_labels(id),''location'',''northeastoutside'')']);
    xlabel('Number of Clusters');
    ylabel('ARI*CalinskiHarabasz index');
    set(gca,'box','on','linewidth',2);
    set(gcf,'Unit','normalized','Position',[0.1,0.1,0.9,0.9]);
    saveas(gcf,[workdir,filesep,'ARI-CH_plot_',dataname{d},'.fig']);
    saveas(gcf,[workdir,filesep,'ARI-CH_plot_',dataname{d},'.png']);
    clear hCH*
    close(figure(1))
    % plot Jaccard_similarity_coef
    figure(1);
    hold on
    
    [f1,xi]=ksdensity(Bestmodel_labels.Jaccard_similarity_coef(:,1)*100,[0:1:100],'BandWidth',2);
    [f2,xi]=ksdensity(Bestmodel_labels.Jaccard_similarity_coef(:,2)*100,[0:1:100],'BandWidth',2);
    fill([xi fliplr(xi)],[zeros(size(f2)),fliplr(f2)],'r',[xi fliplr(xi)],[zeros(size(f1)),fliplr(f1)],'c','facealpha',0.8,'edgecolor','none');
   % plot(xi,f1,'c-',xi,f2,'r-','LineWidth',1);
    legend('Unmatched','Matched','location','northeastoutside');
    xlabel('Jaccard similarity coefficient between shuffles');
    ylabel('PDF');
    set(gca,'box','on','linewidth',2);
    set(gcf,'Unit','normalized','Position',[0.1,0.1,0.9,0.9]);
    saveas(gcf,[workdir,filesep,'Jaccard_similarity_coef_plot_',dataname{d},'.fig']);
    saveas(gcf,[workdir,filesep,'Jaccard_similarity_coef_plot_',dataname{d},'.png']);
    close(figure(1))
    % plot max probability of ensembled labels
    figure(1);
    hold on
    [f1,xi]=ksdensity(Bestmodel_labels.ensembled_label_maxprob(:,1)*100,[50:1:100],'BandWidth',2);
    fill([xi fliplr(xi)],[zeros(size(f1)),fliplr(f1)],'r','facealpha',0.6,'edgecolor','none');
    xlabel('Max probability for ensembled labels');
    ylabel('PDF');
    legend('Max probability','location','northeastoutside');
    set(gca,'box','on','linewidth',2);
    set(gcf,'Unit','normalized','Position',[0.1,0.1,0.9,0.9]);
    saveas(gcf,[workdir,filesep,'max_probability_for_ensembled_labels_',dataname{d},'.fig']);
    saveas(gcf,[workdir,filesep,'max_probability_for_ensembled_labels_',dataname{d},'.png']);
    close(figure(1))
end