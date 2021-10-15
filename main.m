%% CIF AIPW
% ipw estimator
% aipw estimator
% full estimator
% complete case estimator

clear all

% a simulated data from the multiplicative CIF model in Section 5.2
data = csvread('simdata.csv');

%% INPUT
x=data(:,1);
z2=data(:,2); % subject to missingness
z1=data(:,3);
A=data(:,4);
time=data(:,5);
cause=data(:,6);
rselec=data(:,7);
n=size(data,1);

Auxind=1; % with auxiliary variable A

X=[ones(n,1),x];
Z=[z1,z2];
px=size(X,2); % dimension for X including the baseline
pz=size(Z,2); % dimension for Z

% grid points for t
dt=0.1;
times=0.05:dt:3;
ngrid=length(times);

out_case=data(cause==1,:);
out_ntype=data(cause==2,:);
out_censor=data(cause==0,:);

n1=size(out_case,1);
n2=size(out_ntype,1);
n3=size(out_censor,1);

%theta1hat=glmfit(out_case(:,[1,2,5]),out_case(:,7),'binomial');
theta2hat=glmfit(out_ntype(:,[1,2,5]),out_ntype(:,7),'binomial');
ptheta2=length(theta2hat);
theta3hat=glmfit(out_censor(:,[1,2]),out_censor(:,7),'binomial');
ptheta3=length(theta3hat);
ptheta=ptheta2+ptheta3;

%% inverse probability weight
shat=zeros(n,1); % estimated selection prob
dpi=zeros(ptheta,n);
Itheta=zeros(ptheta);
zetai=zeros(ptheta,n);
for i=1:n
    if cause(i)==1
        shat(i)=1;
    elseif cause(i)==2
        Omegai=[1,x(i),z1(i),time(i)];
        shat(i)=expit(Omegai*theta2hat);
        dpi(1:ptheta2,i)=exp(Omegai*theta2hat)...
            /(1+exp(Omegai*theta2hat))^2*Omegai';
        Itheta(1:ptheta2,1:ptheta2)=Itheta(1:ptheta2,1:ptheta2)...
            +exp(Omegai*theta2hat)/(1+exp(Omegai*theta2hat))^2....
            *(Omegai'*Omegai)/n2;
        zetai(1:ptheta2,i)=(rselec(i)-exp(Omegai*theta2hat)...
            /(1+exp(Omegai*theta2hat)))*Omegai'*n/n2;
    else
        Omegai=[1,x(i),z1(i)];
        shat(i)=expit(Omegai*theta3hat);
        dpi((ptheta2+1):ptheta,i)=exp(Omegai*theta3hat)...
            /(1+exp(Omegai*theta3hat))^2*Omegai';
        Itheta((ptheta2+1):ptheta,(ptheta2+1):ptheta)=Itheta((ptheta2+1):ptheta,(ptheta2+1):ptheta)...
            +exp(Omegai*theta3hat)/(1+exp(Omegai*theta3hat))^2....
            *(Omegai'*Omegai)/n3;
        zetai((ptheta2+1):ptheta,i)=(rselec(i)-exp(Omegai*theta3hat)...
            /(1+exp(Omegai*theta3hat)))*Omegai'*n/n3;
    end
end
w_est=rselec./shat;



%% weighted response R
% estimate the censoring distribution by KM
cens=(cause~=0);
[s,tt]=ecdf(time,'Censoring',cens,'function','survivor');
s(end)=[];
tt(1)=[];
Gcx=interp1(tt,s,time,'nearest','extrap');
Gctimes=interp1(tt,s,times,'nearest','extrap');

%% aipw

% initialization
etahatall0=ones(ngrid,px)*0;
gamhat0=ones(1,pz)*0;

[etahatall_a,sigmaeta_a,gamhat_a,sigmagam_a]=esta(n,px,pz,...
    etahatall0,gamhat0,times,ngrid,dt,Gctimes,...
    time,X,Z,A,Auxind,Gcx,cause,w_est,...
    rselec,shat,dpi,Itheta,zetai,ptheta);
%xpd,z1pd1,z1pd2,z2pd1,z2pd2,z2pd3);


%% ipw

% initialization
etahatall0=ones(ngrid,px)*0;
gamhat0=ones(1,pz)*0;

[etahatall_i,sigmaeta_i,gamhat_i,sigmagam_i]=esti(n,px,pz,...
    etahatall0,gamhat0,times,ngrid,dt,Gctimes,...
    time,X,Z,Gcx,cause,w_est,...
    rselec,shat,dpi,Itheta,zetai,ptheta);

save simres.mat

