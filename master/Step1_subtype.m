clc,clear;
data=readtable('G:\brainwander_gitea\jobs\2021_caichao_clustering\PANSS_combat\SZ_GMV_PCA95_data.csv');
outputdir='G:\brainwander_gitea\jobs\2021_caichao_clustering\PANSS_combat';
scriptdir='G:\brainwander_gitea\jobs\2021_caichao_clustering\PANSS_combat';
criteria='CalinskiHarabasz';% 'CalinskiHarabasz','silhouette','DaviesBouldin'
kmeans_out=[outputdir,'\20210923_kmeans_subtype_combat_5fold_',criteria];
hierarchical_out=[outputdir,'\20210923_hierarchical_subtype_combat_5fold_',criteria];
addpath(scriptdir);
subid=data.subid(1:314);
PANSS=table2array(data(1:314,2:31));

%% Optimize Kmean Clustering
if ~exist([kmeans_out,'.mat'],'file')
    XK0=PANSS;
    kparams.subj_num=size(PANSS,1);
    kparams.random_time=10;
    kparams.kmin=2;
    kparams.kstep=1;
    kparams.kmax=10;
    kparams.cvfold=5;
    kparams.kdist={'sqeuclidean','correlation','cityblock','cosine'};
    kparams.kstart={'cluster','uniform','plus'};
    kparams.replicates=100;
    clusters=kparams.kmin:kparams.kstep:kparams.kmax;
    label_kmeans=zeros(kparams.subj_num,length(clusters),length(kparams.kdist),length(kparams.kstart),kparams.random_time,kparams.cvfold);
    CH_kmeans_all=zeros(length(clusters),length(kparams.kdist),length(kparams.kstart),kparams.random_time,kparams.cvfold);
    for c=1:length(clusters)
        for d=1:length(kparams.kdist)
            dist=kparams.kdist{d};
            for s=1:length(kparams.kstart)
                start=kparams.kstart{s};
                for r=1:kparams.random_time
                    part=make_xval_partition(size(XK0,1),kparams.cvfold);
                    for f=1:kparams.cvfold
                        Data=XK0(part~=f,:);
                        kmeans_model=kmeans(Data,clusters(c),'Distance',dist,'Start',start, 'MaxIter',1000,'Replicates',kparams.replicates);
                        label_kmeans(part~=f,c,d,s,r,f)=kmeans_model;
                        eva_kmeans = evalclusters(Data,kmeans_model,criteria);
                        CH_kmeans_all(c,d,s,r,f)=eva_kmeans.CriterionValues;
                    end
                end
                fprintf('cluster #%d  -- dist:%s -- start:%s\n',clusters(c),dist,start);
            end
            
        end
    end
    CH_kmeans=reshape(CH_kmeans_all,length(clusters),length(kparams.kdist),length(kparams.kstart),kparams.random_time*kparams.cvfold);
    CH_kmeans=mean(CH_kmeans,4);
    ARI_kmeans = zeros(length(clusters),length(kparams.kdist),length(kparams.kstart));
    labels=reshape(label_kmeans,[kparams.subj_num,length(clusters),length(kparams.kdist),length(kparams.kstart),kparams.random_time*kparams.cvfold]);
    for c=1:length(clusters)
        for d=1:length(kparams.kdist)
            for s=1:length(kparams.kstart)
                YK_kmeans=squeeze(label_kmeans(:,c,d,s,:));
                tmp_kmeans=cv_cluster_stability(YK_kmeans);
                ARI_kmeans(c,d,s)=tmp_kmeans(1);
            end
        end
    end
    
    [maxCH,I]=max(CH_kmeans(:));
    [CHkmean_c,CHkmean_d,CHkmean_s]=ind2sub(size(CH_kmeans),I);
    subplot(2,2,1);
    plot(clusters,reshape(CH_kmeans,length(clusters),length(kparams.kdist)*length(kparams.kstart)));
    title('CalinskiHarabasz plot');
    ylabel('CalinskiHarabasz value');
    xlabel('number of subtypes');
    
    [maxARI,I]=max(ARI_kmeans(:));
    [ARIkmean_c,ARIkmean_d,ARIkmean_s]=ind2sub(size(ARI_kmeans),I);
    subplot(2,2,2);
    plot(clusters,reshape(ARI_kmeans,length(clusters),length(kparams.kdist)*length(kparams.kstart)));
    title('ARI plot');
    ylabel('ARI value');
    xlabel('number of subtypes');
    
    ARI_CH_kmeans=ARI_kmeans.*CH_kmeans;
    [v,I]=max(ARI_CH_kmeans(:));
    [kmean_c,kmean_d,kmean_s]=ind2sub(size(ARI_CH_kmeans),I);
    subplot(2,2,3);
    plot(clusters,reshape(ARI_CH_kmeans,length(clusters),length(kparams.kdist)*length(kparams.kstart)));    
    title('CalinskiHarabasz*ARI plot');
    ylabel('CalinskiHarabasz*ARI value');
    xlabel('number of subtypes');
    % best model
    Bestmodel={clusters(kmean_c),kparams.kdist{kmean_d},kparams.kstart{kmean_s},ARI_CH_kmeans(kmean_c,kmean_d,kmean_s)};
    best_labels=squeeze(label_kmeans(:,kmean_c,kmean_d,kmean_s,:,:));
    dim=size(best_labels);
    best_labels=reshape(best_labels,dim(1),dim(2)*dim(3));
    Bestmodel_labels=cluster_label_match(best_labels);
    save([kmeans_out,'.mat'],'data','CH_kmeans_all','CH_kmeans','ARI_kmeans','label_kmeans','PANSS','clusters','ARI_CH_kmeans','kparams','Bestmodel','Bestmodel_labels');
    saveas(gcf,[kmeans_out,'_plot.fig']);
end
%% Optimize Hierarchical Clustering
if ~exist([hierarchical_out,'.mat'],'file')
    XK0=PANSS;
    hparams.subj_num=size(PANSS,1);
    hparams.random_time=10;
    hparams.kmin=2;
    hparams.kstep=1;
    hparams.kmax=10;
    hparams.cvfold=5;
    clusters=hparams.kmin:hparams.kstep:hparams.kmax;
    hparams.hdist={'euclidean','squaredeuclidean','seuclidean','mahalanobis','minkowski','chebychev','jaccard','spearman','correlation','cityblock','cosine'};
    hparams.hmethod={'average','centroid','complete','median','single','ward','weighted'};
    label_hier=zeros(hparams.subj_num,length(clusters),length(hparams.hdist),length(hparams.hmethod),hparams.random_time,hparams.cvfold);
    CH_hier_all=zeros(length(clusters),length(hparams.hdist),length(hparams.hmethod),hparams.random_time,hparams.cvfold);
    % for c=1:length(clusters)
    for d=1:length(hparams.hdist)
        dist=hparams.hdist{d};
        for s=1:length(hparams.hmethod)
            method=hparams.hmethod{s};
            for r=1:hparams.random_time
                part=make_xval_partition(size(XK0,1),hparams.cvfold);
                for f=1:hparams.cvfold
                    Data=XK0(part~=f,:);
                    % hierarchical clustering model
                    X=pdist(Data,dist);
                    Z=linkage(X, method);
                    hier_model=cluster(Z,'maxclust',clusters);
                    label_hier(part~=f,:,d,s,r,f)=hier_model;
                    eva_hier =evalclusters(Data,hier_model,criteria);
                    CH_hier_all(:,d,s,r,f)=eva_hier.CriterionValues;
                end
            end
            fprintf('dist:%s -- method:%s\n',dist,method);
        end
        
    end
    % end
    CH_hier=reshape(CH_hier_all,length(clusters),length(hparams.hdist),length(hparams.hmethod),hparams.random_time*hparams.cvfold);
    CH_hier=mean(CH_hier,4);
    ARI_hier = zeros(length(clusters),length(hparams.hdist),length(hparams.hmethod));
    %labels=reshape(label_hier,[hparams.subj_num,length(clusters),length(hparams.hdist),length(hparams.hmethod),hparams.random_time*hparams.cvfold]);
    for c=1:length(clusters)
        for d=1:length(hparams.hdist)
            for s=1:length(hparams.hmethod)
                YK_hier=squeeze(label_hier(:,c,d,s,:));
                tmp_hier=cv_cluster_stability(YK_hier);
                ARI_hier(c,d,s)=tmp_hier(1);
            end
        end
    end
    
    [maxCH,I]=max(CH_hier(:));
    [CHhier_c,CHhier_d,CHhier_s]=ind2sub(size(CH_hier),I);
    subplot(2,2,1);
    plot(clusters,reshape(CH_hier,length(clusters),length(hparams.hdist)*length(hparams.hmethod)));
    title('CalinskiHarabasz plot');
    ylabel('CalinskiHarabasz value');
    xlabel('number of subtypes');
    
    [maxARI,I]=max(ARI_hier(:));
    [ARIhier_c,ARIhier_d,ARIhier_s]=ind2sub(size(ARI_hier),I);
    subplot(2,2,2);
    plot(clusters,reshape(ARI_hier,length(clusters),length(hparams.hdist)*length(hparams.hmethod)));
    title('ARI plot');
    ylabel('ARI value');
    xlabel('number of subtypes');
    
    ARI_CH_hier=ARI_hier.*CH_hier;
    [maxARICH,I]=max(ARI_CH_hier(:));
    [hier_c,hier_d,hier_s]=ind2sub(size(ARI_CH_hier),I);
    subplot(2,2,3);
    plot(clusters,reshape(ARI_CH_hier,length(clusters),length(hparams.hdist)*length(hparams.hmethod)));
    Bestmodel={clusters(hier_c),hparams.hdist{hier_d},hparams.hmethod{hier_s},ARI_CH_hier(hier_c,hier_d,hier_s)};
    title('CalinskiHarabasz*ARI plot');
    ylabel('CalinskiHarabasz*ARI value');
    xlabel('number of subtypes');
    % assemble randsample and fold labels of Bestmodel into one predict
    best_labels=squeeze(label_hier(:,hier_c,hier_d,hier_s,:,:));
    dim=size(best_labels);
    best_labels=reshape(best_labels,dim(1),dim(2)*dim(3));
    Bestmodel_labels=cluster_label_match(best_labels);
    save([hierarchical_out,'.mat'],'data','CH_hier_all','CH_hier','ARI_hier','label_hier','PANSS','clusters','Bestmodel','ARI_CH_hier','hparams','Bestmodel_labels');
    saveas(gcf,[hierarchical_out,'_plot.fig']);
end




