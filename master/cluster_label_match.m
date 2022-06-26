function output=cluster_label_match(umatched_labels,nperm)
% umatched_labels:[nsub nfolds] labels
% lab_merge:merged label
% load hierarchical_subtype_results;
% umatched_labels=best_label_hier_nfold;
m=size(umatched_labels);
matched_labels=zeros(m);
ref=umatched_labels(:,1);
matched_labels(:,1)=ref;
Jaccard_similarity_coef=zeros(m(2),1);
labnum=unique(umatched_labels);
labnum(labnum==0)=[];
%sort label
dice_orig=[];
dice_new=[];
for k=1:m(2)
    targ=umatched_labels(:,k);
    nonzero_idx=all([ref targ],2);
    ref1=ref(nonzero_idx);
    targ1=targ(nonzero_idx);
    Jaccard_similarity_coef(k,1)=1-pdist([ref1,targ1]','jaccard');
    dice_orig(:,:,k)=dice_matrix(ref1,targ1,labnum);    
    perm_matrix=[];
    if length(labnum)<=6  % exact permutation
        perm_matrix=perms(labnum);
    else  % randperm
        for k=1:nperm
            idx=randperm(length(labnum));
            perm_matrix(k,:)=labnum(idx);
        end
    end
    mean_maxdice=0;
    newtarg=[];
    for iter=1:size(perm_matrix,1)
        tmplab=perm_matrix(iter,:);
        tfid=[];
        maxdice=[]; 
        tmptarg=zeros(size(targ));
        for i=1:length(tmplab)
            rlab=tmplab(i);
            tmpdice=dice_orig(rlab,:,k);
            [sortdice,sid]=sort(tmpdice,'descend');   
            idx=find(~ismember(sid,tfid));
            tfid(i)=sid(idx(1));
            maxdice(i)=sortdice(idx(1));  
            tmptarg(targ==tfid(i))=rlab ;
         end
        if mean(maxdice)>mean_maxdice
           newtarg=tmptarg;
           mean_maxdice=mean(maxdice); 
        end 
        if mod(iter,500)==1
           fprintf('fold #%d - iter #%d - mean dice:%0.6f\n',k,iter,mean_maxdice);
        end
    end  
    dice_new(:,:,k)=dice_matrix(ref1,newtarg(nonzero_idx),labnum);    
    matched_labels(:,k)=newtarg;
    Jaccard_similarity_coef(k,2)=1-pdist([ref1,newtarg(nonzero_idx)]','jaccard');
end

%calcualte the label probability for each subject
lab_prob=zeros(m(1),length(labnum));
for k=1:m(1)
    tmp=matched_labels(k,:);
    t=sum(tmp>0);
    for n=1:length(labnum)
        nlab=sum(tmp==labnum(n));
        lab_prob(k,n)=nlab/t;
    end
end
[maxprob,maxprob_label]=max(lab_prob,[],2);
output.matched_labels=matched_labels;
output.umatched_labels=umatched_labels;
output.ensembled_label=maxprob_label;
output.ensembled_label_maxprob=maxprob;
output.Jaccard_similarity_coef=Jaccard_similarity_coef;
output.dice_matrix_new=dice_new;
output.dice_matrix_orig=dice_orig;
end

%% subfunctions 
function M=dice_matrix(ref,targ,labnum)
    for i=1:length(labnum)
        row=ref==labnum(i);
        for j=1:length(labnum)
          col=targ==labnum(j);
          conj=row.*col;
          M(i,j)=2*sum(conj)/(sum(row)+sum(col));
        end
    end
end

% function matched_label=clustering_label_match();