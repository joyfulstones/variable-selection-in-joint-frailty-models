%macro nlm(data);
proc nlmixed data=&data qpoints=10;
parms beta1=1  beta6=-1 gamma=1 alpha1=1 
alpha2=-0.5  theta=1 
h1=5 h2=5 h3=5 h4=5 h5=5 h6=5 h7=5 h8=5 h9=5 h10=5 
r1=8 r2=8 r3=8 r4=8 r5=8 r6=8 r7=8 r8=8 r9=8 r10=8; 

bounds theta >0, h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10>=0;
if Ddelta=0 then loglik=beta1*Z1+beta6*Z6
+nu+log(pr1*r1+pr2*r2+pr3*r3+pr4*r4+pr5*r5+pr6*r6+pr7*r7+pr8*r8+pr9*r9+pr10*r10);
if Ddelta=1 then do;
temp1=-exp(beta1*Z1+beta6*Z6+nu)
*(r1*Lr1+r2*Lr2+r3*Lr3+r4*Lr4+r5*Lr5+r6*Lr6+r7*Lr7+r8*Lr8+r9*Lr9+r10*Lr10);
if delta=0 then temp2=0;
if delta=1 then temp2=alpha1*Z1+alpha2*Z2
+gamma*nu+log(p1*h1+p2*h2+p3*h3+p4*h4+p5*h5+p6*h6+p7*h7+p8*h8+p9*h9+p10*h10);
temp3=-exp(alpha1*Z1+alpha2*Z2+gamma*nu)
*(h1*L1+h2*L2+h3*L3+h4*L4+h5*L5+h6*L6+h7*L7+h8*L8+h9*L9+h10*L10);
temp4=1/theta*(nu-exp(nu))-1/theta*log(theta)-lgamma(1/theta);
temp5=nu*nu/2;
loglik=temp1+temp2+temp3+temp4+temp5;
end;

model X~ general(loglik);
random nu ~normal([0],[1]) subject=id;
ods output ParameterEstimates=two.test1;
run;
%mend();
%macro estimate();
libname two 'C:\Users\28029\Desktop\code\SAS\temp  result';
%do jj=1%to 1;
proc import out=two.recurrent_event&jj 
datafile="C:\Users\28029\Desktop\code\SAS\data\recurrent_event&jj..csv" 
dbms=csv
replace;
run;

data two.one;
set two.recurrent_event&jj;
if delta=1 and Ddelta=1;
run;

proc capability data=two.one noprint;
var X; 
output out=two.quant_death pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
run;

data two.quant_death;
set two.quant_death;
aa=1;
run;

data two.recurrent_event&jj;
set two.recurrent_event&jj;
aa=1;
run;

data two.recurrent_event&jj;
merge two.recurrent_event&jj two.quant_death;
by aa;
run;


data two.oner;
set two.recurrent_event&jj;
if  Ddelta=0;
run;

proc capability data=two.oner noprint;
var T; 
output out=two.quant_recurrent pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr; 
run;

data two.quant_recurrent;
set two.quant_recurrent;
bb=1;
run;

data two.recurrent_event&jj;
set two.recurrent_event&jj;
bb=1;
run;

data two.recurrent_event&jj;
merge two.recurrent_event&jj two.quant_recurrent;
by bb;
run;







data two.recurrent_event&jj;
set two.recurrent_event&jj;
array p {10} p1-p10;
array L {10} L1-L10;
array pr {10} pr1-pr10;
array Lr {10} Lr1-Lr10;
array q{11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;
array tr{11} qr0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;


do ii=1 to 10;
	p{ii}=0;
	L{ii}=0;
    pr{ii}=0;
	Lr{ii}=0;
end;

if X>=q{1} then do;
p{1}=1;
end;
do ii=2 to 10;
if X>=q{ii} then do;
p{ii-1}=0;
p{ii}=1;
end;
end;
if X>q{11} then do;
p{10}=0;
end;

if T>=tr{1} then do;
pr{1}=1;
end;
do ii=2 to 10;
if T>=tr{ii} then do;
pr{ii-1}=0;
pr{ii}=1;
end;
end;
if T>tr{11} then do;
pr{10}=0;
end;


do ii=1 to 10;
if X>=q{ii} then do;
if X>=q{ii+1} then do;
L{ii}=q{ii+1}-q{ii};
end;
else do;
L{ii}=X-q{ii};
end;
end;
end;


do ii=1 to 10;
if X>=tr{ii} then do;
if X>=tr{ii+1} then do;
Lr{ii}=tr{ii+1}-tr{ii};
end;
else do;
Lr{ii}=X-tr{ii};
end;
end;
end;
run;



%nlm(two.recurrent_event&jj) 

data two.one&jj;
set _null_;
run;
data two.one&jj;
set two.one&jj two.test1;
KEEP parameter Estimate;
run;

data two.quant_death&jj;
set two.quant_death(drop=aa);
run;

data two.quant_recurrent&jj;
set two.quant_recurrent(drop=bb);
run;
%end;
%mend();
%estimate()

%macro toge();
data two.estall;  
set _null_;
run;

data two.quant_deathall;
set _null_;
run;

data two.quant_recurrentall;
set _null_;
run;


%do ii=1%to 1;
data two.estall;   
set two.estall two.one&ii;   
run;
data two.quant_deathall;
set two.quant_deathall two.quant_death&ii;
run;

data two.quant_recurrentall;
set two.quant_recurrentall two.quant_recurrent&ii;
run;
%end;

proc export data=two.estall
outfile="C:\Users\28029\Desktop\code\SAS\result\estall.csv"
dbms=csv
label
replace;
run;

proc export data=two.quant_deathall
outfile="C:\Users\28029\Desktop\code\SAS\result\quant_deathall.csv"
dbms=csv
label
replace;
run;

proc export data=two.quant_recurrentall
outfile="C:\Users\28029\Desktop\code\SAS\result\quant_recurrentall.csv"
dbms=csv
label
replace;
run;
%mend();
%toge()
