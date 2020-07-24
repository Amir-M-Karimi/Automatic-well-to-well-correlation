function [ownew] = PCA_observation_attributes(Avg,Cv,R,theta);
numdata=length(Avg); %should be even
%step 1, generating a dataset

x=Avg;
y=Cv;
%step 2, finding a mean and subtracting
xmean=mean(x);
ymean=mean(y);
xnew=x-xmean*ones(numdata,1);
ynew=y-ymean*ones(numdata,1);
covariancematrix=cov(xnew,ynew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));
finaldata1=maxeigval'*[xnew,ynew]';
ownew1=finaldata1;

z=R;
s=theta;
%step 2, finding a mean and subtracting
zmean=mean(z);
smean=mean(s);
snew=s-smean*ones(numdata,1);
znew=z-zmean*ones(numdata,1);
covariancematrix=cov(znew,snew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));
finaldata2=maxeigval'*[znew,snew]';
ownew2=finaldata2;

x=ownew1';
y=ownew2';
%step 2, finding a mean and subtracting
xmean=mean(x);
ymean=mean(y);
xnew=x-xmean*ones(numdata,1);
ynew=y-ymean*ones(numdata,1);
covariancematrix=cov(xnew,ynew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));
finaldata=maxeigval'*[xnew,ynew]';
ownew=finaldata;
end

