function [score,stdscore]=cv_cluster_stability(S)
k=0;
for i=1:size(S,2)-1
    for j=i+1:size(S,2)
        k=k+1;
        zero_idx=any([S(:,i) S(:,j)]==0,2);
        [a(k),b(k),c(k),d(k)]=RandIndex(S(~zero_idx,i),S(~zero_idx,j));
    end
end
score=[mean(a) mean(b) mean(c) mean(d)];
stdscore=[std(a) std(b) std(c) std(d)];
end