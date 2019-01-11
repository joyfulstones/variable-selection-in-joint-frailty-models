function betav=fbeta(Z,X,delta,Lambdav,Lambda_posi,gamma0v,gamma1v,mindeathX,STK)
%{
temp=sum(Z.*delta+log(exp(-gamma0v-gamma1v.*Z)./(exp(-gamma0v-gamma1v.*Z)+1)).*Z);
f=@(x)temp-sum(log(exp(Lambdav(Lambda_posi)+x*Z.*mindeathX-gamma0v-gamma1v.*Z)./(exp(Lambdav(Lambda_posi)+x*Z.*mindeathX-gamma0v-gamma1v.*Z)+1)).*Z);
opt=optimset('Display','off');
betav=fsolve(f,0.5,opt);
%}

%{
temp=sum(Z.*delta+log(exp(-gamma0v-gamma1v.*Z)./(exp(-gamma0v-gamma1v.*Z)+1)).*Z);
f=@(x)temp-sum(log(exp(Lambdav(Lambda_posi)+x*Z.*X-gamma0v-gamma1v.*Z)./(exp(Lambdav(Lambda_posi)+x*Z.*X-gamma0v-gamma1v.*Z)+1)).*Z);
opt=optimset('Display','off');
betav=fsolve(f,0.5,opt);
%}

temp=sum(Z.*delta+log(exp(-gamma0v-gamma1v.*Z)./(exp(-gamma0v-gamma1v.*Z)+1)).*Z);
f=@(x)temp-sum(STK.*log(exp(Lambdav(Lambda_posi)+x*Z.*X-gamma0v-gamma1v.*Z)./(exp(Lambdav(Lambda_posi)+x*Z.*X-gamma0v-gamma1v.*Z)+1)).*Z);
opt=optimset('Display','off');
betav=fsolve(f,0.5,opt);
