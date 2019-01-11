function Lambdav=fLambda(Z,death_time,death_number,ST,S2T,betav,gamma0v,gamma1v,X)
Lambdav=zeros(death_number+1,1);
temp1=exp(-gamma0v-gamma1v.*Z);
temp=1+sum(ST(:,1).*log(temp1./(1+temp1)));
temp2=sum(S2T(:,1).*log(temp1./(1+temp1)));
f=@(x)sum(ST(:,1).*log((exp(x+betav.*Z.*death_time(1)).*temp1)./(1+exp(x+betav.*Z.*death_time(1)).*temp1)))+...
sum(S2T(:,1).*log((exp(betav.*Z.*X).*temp1)./(1+exp(betav.*Z.*X).*temp1)))-temp-temp2;
opt=optimset('Display','off');
Lambdav(2)=fsolve(f,0,opt);
for i=3:death_number+1
    temp=1+sum(ST(:,i-1).*log(temp1./(exp(-Lambdav(i-1)-betav.*Z.*death_time(i-2))+temp1)));
    temp2=sum(S2T(:,i-1).*log(temp1./(exp(-Lambdav(i-1)-betav.*Z.*death_time(i-2))+temp1)));
    f=@(x)sum(ST(:,i-1).*log(temp1./(exp(-x-betav.*Z.*death_time(i-1))+temp1)))+...
    sum(S2T(:,i-1).*log(temp1./(exp(-Lambdav(i-1)-betav.*Z.*X)+temp1)))-temp-temp2;
    Lambdav(i)=fsolve(f,0,opt);
end

%{
Lambdav=zeros(death_number+1,1);
temp1=exp(-gamma0v-gamma1v.*Z);
temp=1+sum(ST(:,1).*log((exp(betav.*Z.*death_time(1)).*temp1)./(1+exp(betav.*Z.*death_time(1)).*temp1)));
f=@(x)sum(ST(:,1).*log((exp(x+betav.*Z.*death_time(1)).*temp1)./(1+exp(x+betav.*Z.*death_time(1)).*temp1)))-temp;
opt=optimset('Display','off');
Lambdav(2)=fsolve(f,0,opt);
for i=3:death_number+1
    temp=1+sum(ST(:,i-1).*log(temp1./(exp(-Lambdav(i-1)-betav.*Z.*death_time(i-1))+temp1)));
    f=@(x)sum(ST(:,i-1).*log(temp1./(exp(-x-betav.*Z.*death_time(i-1))+temp1)))-temp;
    Lambdav(i)=fsolve(f,0,opt);
end
%}


