tic;
s=RandStream('mt19937ar','seed','shuffle'); 
RandStream.setGlobalStream(s);
n=500;
B=1;
for i=1:B
   temp1=strcat('recurrent_event',num2str(i),'.csv');
num1=csvread(temp1,1,0);
for j=1:1
    id=num1(:,9);
   temp_data=zeros(10000,12);
   temp_id=unidrnd(n,n,1);
   count_row=0;
   temp_begin=1;
   temp2_id=zeros(10000,1);
   for K=1:n
       temp_number=find(id==temp_id(K));
       count_row=count_row+sum(temp_number~=0);
       temp_data(temp_begin:count_row,:)=num1(temp_number,:);
       temp2_id(temp_begin:count_row)=K;
       temp_begin=count_row+1;       
   end
   temp_data=temp_data(1:count_row,:);
   Z1=temp_data(:,1);
   Z2=temp_data(:,2);
   Z3=temp_data(:,3);
   Z4=temp_data(:,4);
   Z5=temp_data(:,5);
   Z6=temp_data(:,6);
   X=temp_data(:,7);
   T=temp_data(:,8);
   delta=temp_data(:,10);
   Ddelta=temp_data(:,11);
   N=max(size(Ddelta))-n;
   N=repmat(N,count_row,1);
   id=temp2_id(temp2_id~=0);
columns ={'Z1','Z2','Z3','Z4','Z5','Z6','X','T','id','delta','Ddelta','N'};
data=table(Z1,Z2,Z3,Z4,Z5,Z6,X,T,id,delta,Ddelta,N,'VariableNames',columns);
filename=sprintf('recurrent_event%db%d.csv',i,j);
writetable(data,filename);
id=num1(:,9); %idҪ���ԭ����
end  
end
toc;
