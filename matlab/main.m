s=RandStream('mt19937ar','seed','shuffle'); 
RandStream.setGlobalStream(s);
n=500;
alpha=[1.0,-0.5,0,0,0,0]';
beta=[1.0,0,0,0,0,-1.0]';
gamma=1.0;
theta=1.0;
B=1;
for i=1:B
    [Z,X,T,id,delta,Ddelta,N]=tempdata(alpha,beta,gamma,theta,n);
    Z1=Z(:,1);
    Z2=Z(:,2);
    Z3=Z(:,3);
    Z4=Z(:,4);
    Z5=Z(:,5);
    Z6=Z(:,6);
    N=repmat(N,N+n,1);

    
columns ={'Z1','Z2','Z3','Z4','Z5','Z6','X','T','id','delta','Ddelta','N'};
data1=table(Z1,Z2,Z3,Z4,Z5,Z6,X,T,id,delta,Ddelta,N,'VariableNames', columns);
filename=sprintf('recurrent_event%d.csv',i);
writetable(data1,filename);

i
end
