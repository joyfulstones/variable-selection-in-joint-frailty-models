tic;
s=RandStream('mt19937ar','seed','shuffle'); 
RandStream.setGlobalStream(s);
warning('off');
n=400;
beta=0.5;
gamma0=1;
gamma1=-0.5;
B=1000;
betaest=zeros(B,1);
gamma0est=zeros(B,1);
gamma1est=zeros(B,1);
Lambdaest=zeros(B,n);
SEE_beta=zeros(B,1);
SEE_gamma0=zeros(B,1);
SEE_gamma1=zeros(B,1);
det=zeros(B,1);
A=zeros(3*B,3);
ZL=zeros(B,n);
deltaL=zeros(B,n);
XL=zeros(B,n);
i=1;
K=0;
LL=0;
while (i<=B)
    re=est(beta,gamma0,gamma1,n);
    betaest(i)=re(1);
    gamma0est(i)=re(2);
    gamma1est(i)=re(3);
    death_number=re(4);
    Z=re(death_number+6:death_number+5+n);
    X=re(death_number+6+n:death_number+5+n*2);
    delta=re(death_number+6+n*2:death_number+5+n*3);
    Lambda_posi=re(death_number+6+n*3:death_number+5+n*4);
    death_time=re(death_number+6+n*4:death_number+5+n*4+sum(delta));
    Lambdaest(i,1:sum(delta)+1)=(re(5:death_number+5))';
    [SEE_beta(i),SEE_gamma0(i),SEE_gamma1(i),A((i-1)*3+1:i*3,:)]=vari(betaest(i),gamma0est(i),gamma1est(i),re(5:death_number+5),Z,X,delta,n,Lambda_posi,death_time,death_number);   
    tempSEE=[SEE_beta(i),SEE_gamma0(i),SEE_gamma1(i)];
    if(max(tempSEE)<1e4)
     i
    i=i+1;
    end
    if(max(tempSEE)>=1e4)
        K=K+1;
    end
    tempLL=sum(isnan(tempSEE));
    if(tempLL~=0)
        LL=LL+1;
    end
end
cp_beta=mean((repmat(beta,B,1)>=(betaest-1.96.*sqrt(SEE_beta./n))).*(repmat(beta,B,1)<=(betaest+1.96.*sqrt(SEE_beta./n))));
cp_gamma0=mean((repmat(gamma0,B,1)>=(gamma0est-1.96.*sqrt(SEE_gamma0./n))).*(repmat(gamma0,B,1)<=(gamma0est+1.96.*sqrt(SEE_gamma0./n))));
cp_gamma1=mean((repmat(gamma1,B,1)>=(gamma1est-1.96.*sqrt(SEE_gamma1./n))).*(repmat(gamma1,B,1)<=(gamma1est+1.96.*sqrt(SEE_gamma1./n))));
se_beta=sqrt(var(betaest));
se_gamma0=sqrt(var(gamma0est));
se_gamma1=sqrt(var(gamma1est));
ASE_beta=mean(sqrt(SEE_beta./n));
ASE_gamma0=mean(sqrt(SEE_gamma0./n));
ASE_gamma1=mean(sqrt(SEE_gamma1./n));
AE_beta=mean(betaest);
AE_gamma0=mean(gamma0est);
AE_gamma1=mean(gamma1est);
toc;