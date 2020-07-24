function [ wwnew,wwnew1,wwnew2] = PCA_witness_logs(ww1,ww2,ww3,ww4);
numdata=length(ww1); %should be even
%step 1, generating a dataset
x=ww1;
y=ww2;
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
wwnew1=finaldata1;

z=ww3;
s=ww4;
%step 2, finding a mean and subtracting
zmean=mean(z);
smean=mean(s);
znew=z-zmean*ones(numdata,1);
snew=s-smean*ones(numdata,1);
covariancematrix=cov(znew,snew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));
finaldata2=maxeigval'*[znew,snew]';
wwnew2=finaldata2;

x=wwnew1';
y=wwnew2';
%step 2, finding a mean and subtracting
xmean=mean(x);
ymean=mean(y);
xnew=x-xmean*ones(numdata,1);
ynew=y-ymean*ones(numdata,1);
covariancematrix=cov(xnew,ynew);
[V,D] = eig(covariancematrix);
D=diag(D);
maxeigval=V(:,find(D==max(D)));
finaldata3=maxeigval'*[xnew,ynew]';
wwnew=finaldata3;








end

