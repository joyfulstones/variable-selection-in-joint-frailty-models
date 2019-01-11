function [SEE_betaest,SEE_gamma0est,SEE_gamma1est,A]=vari(betav,gamma0v,gamma1v,Lambdav,Z,X,delta,n,Lambda_posi,death_time,death_number)
LambdaX=Lambdav(Lambda_posi);%这里大于最后一个死亡点的并没有令为无穷，后面编程要当成无穷
GX=exp(LambdaX+betav.*Z.*X-gamma0v-gamma1v.*Z)./(1+exp(LambdaX+betav.*Z.*X-gamma0v-gamma1v.*Z)).*(X<=death_time(death_number))+(X>death_time(death_number));
WX=[Z';GX';(GX.*Z)'];
HX=[Z';(GX.*(1-delta))';(GX.*Z.*(1-delta))'];
YX=(repmat(X,1,n)>=repmat(X',n,1));
GX2=exp((repmat(LambdaX',n,1))+betav.*repmat(Z,1,n).*repmat(X',n,1)-gamma0v-gamma1v.*repmat(Z,1,n))./(1+exp((repmat(LambdaX',n,1))+betav.*repmat(Z,1,n).*repmat(X',n,1)...
     -gamma0v-gamma1v.*repmat(Z,1,n))).*((repmat(X',n,1))<=death_time(death_number))+((repmat(X',n,1))>death_time(death_number));
BXX=zeros(n,n);%B(X_i,X_j)
temp1=(sum(YX.*(1-GX2)));
temp1=temp1+(temp1==0);
temp=delta.*GX./temp1';
for i=1:n
    temp2=(repmat(temp,1,n)).*(((repmat(X,1,n))<X(i)).*(repmat(X,1,n)>repmat(X',n,1)).*(X(i)>=(repmat(X',n,1)))...
            -((repmat(X,1,n))>X(i)).*(repmat(X,1,n)<repmat(X',n,1)).*(X(i)<(repmat(X',n,1))));
      BXX(i,:)=-(sum(temp2));
    BXX(i,:)=exp(-(sum(temp2)));
end
hX=zeros(3,n);
temp4=(1-repmat(GX,1,n)).*YX;
for i=1:n
    for j=1:n
        if(temp4(i,j)==0)
            BXX(j,i)=0;
        end
    end
end
temp3=(1-repmat(GX,1,n)).*YX.*BXX';
for i=1:3
hX(i,:)=(sum(repmat(HX(i,:)',1,n).*temp3))./temp1;
end
Sigma=(WX-hX)*((WX-hX)'.*repmat(delta,1,3))./n;
A=zeros(3,3);
A(:,[2,3])=(WX-hX)*([GX,Z.*GX].*repmat(delta,1,2))./n;

%{
A(:,1)=(WX-hX)*(-Z.*X.*GX.*delta)./n+(WX-hX)*(Z.*delta./(betav.*Z+2.*X))./n;

varia=inv(A)*Sigma*(inv(A))';
SEE_betaest=varia(1,1);
SEE_gamma0est=varia(2,2);
SEE_gamma1est=varia(3,3);
%}

sort_d=sort(death_time);
qua=0:0.1:1;
quaX=quantile(sort_d,qua);
quaX=quaX';
quaXL=[0;quaX];
dis=quaXL(2:12)-quaXL(1:11);
lambda=zeros(11,1);
Lambda1=0;
for i=1:10
    tempposi=sum(quaX(i)==death_time);
    if(tempposi==0)
        tempposi=sum(quaX(i)>=death_time);
        tempposi=tempposi+1;
        Lambda2=Lambdav(tempposi);
        lambda(i)=(Lambda2-Lambda1)./dis(i);
        Lambda1=Lambda2;
    else
      tempposi=sum(quaX(i)>=death_time);
        Lambda2=Lambdav(tempposi);
        lambda(i)=(Lambda2-Lambda1)./dis(i);
        Lambda1=Lambdav(tempposi+1);
        
    end
end
lambda(11)=(Lambdav(death_number+1)-Lambda1)./dis(11);
lambdaX_posi=(sum((repmat(X',11,1))>=(repmat(quaXL(1:11),1,n))));
lambdaX=lambda(lambdaX_posi);
A(:,1)=(WX-hX)*(-Z.*X.*GX.*delta)./n+(WX-hX)*(Z.*delta./(betav.*Z+lambdaX+(lambdaX==0)))./n;
varia=inv(A)*Sigma*(inv(A))';
SEE_betaest=varia(1,1);
SEE_gamma0est=varia(2,2);
SEE_gamma1est=varia(3,3);



%{

tempinteg=zeros(n,3);
WD=zeros(3,death_number);
hD=zeros(3,death_number);
for i=1:death_number
    WD(:,i)=(sum(WX'.*repmat((death_time(i)==X),1,3)))';
    hD(:,i)=(sum(hX'.*repmat((death_time(i)==X),1,3)))';
end
dis=death_time(2:death_number)-death_time(1:death_number-1);
dis=[death_time(1);dis];
for i=1:n
int_posi=sum(X(i)>=death_time)+1;
if(int_posi==1)
tempinteg(i,:)=((WD(:,1)-hD(:,1)).*X(i))';
elseif(int_posi==(death_number+1))
    tempinteg(i,:)=(sum((WD(:,1:int_posi-1)-hD(:,1:int_posi-1))'.*(repmat(dis,1,3))));
else
    tempinteg(i,:)=((WD(:,int_posi)-hD(:,int_posi)).*(X(i)-death_time(int_posi-1)))'+(sum((WD(:,1:int_posi-1)-hD(:,1:int_posi-1))'.*(repmat(dis(1:int_posi-1),1,3))));
end
end
integ=mean(tempinteg.*repmat(Z,1,3));
A(:,1)=(WX-hX)*(-Z.*X.*GX.*delta)./n+integ';
varia=inv(A)*Sigma*(inv(A))';
det=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
SEE_betaest=varia(1,1);
SEE_gamma0est=varia(2,2);
SEE_gamma1est=varia(3,3);
%}

%{
pointnumber=1001;
sort_X=sort(X);
XT=linspace(0,sort_X(n),pointnumber);%当实验时间较长时可以增加分点的个数，可以考虑增加到1000
XT=XT';
XT=XT(2:pointnumber);
pointnumber=pointnumber-1;
LambdaT_posi=sum(repmat(XT',death_number,1)>=repmat(death_time,1,pointnumber));
LambdaT_posi=LambdaT_posi+1;
LambdaXT=Lambdav(LambdaT_posi);
dis=XT(2)-XT(1);
BTX=zeros(pointnumber,n);%B(t_i,X_j)
for i=1:pointnumber
    temp4=(repmat(temp,1,n)).*(((repmat(X,1,n))<XT(i)).*(repmat(X,1,n)>repmat(X',n,1)).*(XT(i)>=(repmat(X',n,1)))...
            -((repmat(X,1,n))>XT(i)).*(repmat(X,1,n)<repmat(X',n,1)).*(XT(i)<(repmat(X',n,1))));
    BTX(i,:)=exp(-(sum(temp4)));
end

YXT=(repmat(X,1,pointnumber)>=repmat(XT',n,1));
GXT2=exp((repmat(LambdaXT',n,1))+betav.*repmat(Z,1,pointnumber).*repmat(XT',n,1)-gamma0v-gamma1v.*repmat(Z,1,pointnumber))./(1+exp((repmat(LambdaXT',n,1))+betav.*repmat(Z,1,pointnumber).*repmat(XT',n,1)...
     -gamma0v-gamma1v.*repmat(Z,1,pointnumber))).*((repmat(XT',n,1))<=death_time(death_number))+((repmat(XT',n,1))>death_time(death_number));

temp5=(sum(YXT.*(1-GXT2)));
temp5=temp5+(temp5==0);
hXT=zeros(3,pointnumber);
temp6=(1-repmat(GX,1,pointnumber)).*YXT.*BTX';
for i=1:3
hXT(i,:)=(sum(repmat((HX(i,:))',1,pointnumber).*temp6))./temp5;
end

tempinteg=zeros(n,3);
for i=1:n
W=[repmat(Z(i),1,pointnumber);exp(LambdaXT'+betav.*Z(i).*XT'-gamma0v-gamma1v.*Z(i))./(1+exp(LambdaXT'+betav.*Z(i).*XT'-gamma0v-gamma1v.*Z(i)))...
    .*(XT<=death_time(death_number))'+(XT>death_time(death_number))';
    Z(i).*exp(LambdaXT'+betav.*Z(i).*XT'-gamma0v-gamma1v.*Z(i))./(1+exp(LambdaXT'+betav.*Z(i).*XT'-gamma0v-gamma1v.*Z(i)))...
    .*(XT<=death_time(death_number))'+Z(i).*(XT>death_time(death_number))'];
int_posi=sum(X(i)>=XT)+1;
if(int_posi==1)
tempinteg(i,:)=((W(:,1)-hXT(:,1)).*X(i))';
elseif(int_posi==(pointnumber+1))
    tempinteg(i,:)=(sum((W(:,1:int_posi-1)-hXT(:,1:int_posi-1))'.*dis));
else
    tempinteg(i,:)=((W(:,int_posi)-hXT(:,int_posi)).*(X(i)-XT(int_posi-1)))'+(sum((W(:,1:int_posi-1)-hXT(:,1:int_posi-1))'.*dis));
end
end
integ=mean(tempinteg.*repmat(Z,1,3));
A(:,1)=(WX-hX)*(-Z.*X.*GX.*delta)./n+integ';
varia=inv(A)*Sigma*(inv(A))';
SEE_betaest=varia(1,1);
SEE_gamma0est=varia(2,2);
SEE_gamma1est=varia(3,3);
%}


