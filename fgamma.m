function gammav=fgamma(Z,X,delta,Lambdav,Lambda_posi,gamma0v,gamma1v,betav,STK)
temp1=sum((delta+(1-delta)./(1+exp(Lambdav(Lambda_posi)+betav.*Z.*X-gamma0v-gamma1v.*Z))).*STK);
temp1=temp1+sum(delta.*(1-STK));
temp2=sum((Z.*(delta+(1-delta)./(1+exp(Lambdav(Lambda_posi)+betav.*Z.*X-gamma0v-gamma1v.*Z)))).*STK);
temp2=temp2+sum(delta.*(1-STK).*Z);
f=@(x)[temp1-sum(exp(x(1)+x(2)*Z)./(1+exp(x(1)+x(2)*Z)));temp2-sum((exp(x(1)+x(2)*Z)./(1+exp(x(1)+x(2)*Z))).*Z)];
opt=optimset('Display','off');
gammav=fsolve(f,[1,-1],opt);