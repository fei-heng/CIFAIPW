function  [etahatall,sigmaeta,gamhat,sigmagam]=estf(n,px,pz,...
    etahatall0,gamhat0,times,ngrid,dt,Gctimes,...
    time,X,Z,Gcx,cause)

delta=(cause==1);

F1hat=zeros(n,ngrid);
Rhat=zeros(n,ngrid);
Deta=zeros(n,px,ngrid);
Dgam=zeros(n,pz,ngrid);
%Wt=zeros(n,n,ngrid);
%Psihat=zeros(n,n,ngrid);

Beta=zeros(px,ngrid);
Bgam=zeros(pz,ngrid);
BKB=zeros(pz,ngrid);

Cetaeta=zeros(px,px,ngrid);
Cgameta=zeros(pz,px,ngrid);
Cgamgam=zeros(pz,pz,ngrid);
K=zeros(pz,px,ngrid);
CKC=zeros(pz,pz,ngrid);

% estf function
etahatall=etahatall0;
gamhat=gamhat0;
error = 1;
iter = 1;
while error>0.00001 && iter<30
    for igrid=1:ngrid
        t=times(igrid);
        etahat=etahatall(igrid,:);
        
        Rhat(:,igrid)=(time<=t).*delta./Gcx;
        F1hat(:,igrid)=1-exp(-(X*etahat').*exp(Z*gamhat'));
        Deta(:,:,igrid)=repmat(exp(-(X*etahat').*exp(Z*gamhat')).*exp(Z*gamhat'),1,px).*X;
        Dgam(:,:,igrid)=repmat(exp(-(X*etahat').*exp(Z*gamhat')).*(X*etahat').*exp(Z*gamhat'),1,pz).*Z;
        %Wt(:,:,igrid)=eye(n);
        %Psihat(:,:,igrid)=diag(w_est);
        
        Beta(:,igrid)=Deta(:,:,igrid)'*(Rhat(:,igrid)-F1hat(:,igrid));
        Bgam(:,igrid)=Dgam(:,:,igrid)'*(Rhat(:,igrid)-F1hat(:,igrid));
        
        
        Cetaeta(:,:,igrid)=Deta(:,:,igrid)'*Deta(:,:,igrid);
        Cgameta(:,:,igrid)=Dgam(:,:,igrid)'*Deta(:,:,igrid);
        Cgamgam(:,:,igrid)=Dgam(:,:,igrid)'*Dgam(:,:,igrid);
        
        K(:,:,igrid)=Cgameta(:,:,igrid)*pinv(Cetaeta(:,:,igrid));
        BKB(:,igrid)=Bgam(:,igrid)-K(:,:,igrid)*Beta(:,igrid);
        CKC(:,:,igrid)=Cgamgam(:,:,igrid)-K(:,:,igrid)*Cgameta(:,:,igrid)';
    end
    
    Igaminv=pinv(sum(CKC,3)*dt);
    BKBint=sum(BKB,2)*dt;
    gamchange=(Igaminv*BKBint)';
    gamhat1=gamhat+gamchange;
    
    etachange=zeros(ngrid,px);
    for igrid=1:ngrid
        etachange(igrid,:)=(pinv(Cetaeta(:,:,igrid))*(Beta(:,igrid)...
            -(Cgameta(:,:,igrid))'*(gamchange)'))';
    end
    etahatall1=etahatall+etachange;
    
    error=max(abs(gamchange),[],'all')+max(abs(etachange),[],'all');
    etahatall=etahatall1;
    gamhat=gamhat1;
    
    iter = iter + 1;
end

if iter == 30
    etahatall=nan(ngrid,px);
    gamhat=nan(1,pz);
end

%% variance estimation
% Phieta, PhiG
Phietaihat=zeros(px,n,ngrid);
PhiGihat=zeros(pz,n);
qetahat=zeros(px,ngrid,ngrid);
qgamhat=zeros(pz,ngrid,ngrid);
qkq=zeros(pz,ngrid);
logGc=-log(Gctimes);
dlogGc=logGc-[0,logGc(1:(ngrid-1))];

for sgrid=1:ngrid
    for igrid=1:ngrid
        qetahat(:,sgrid,igrid)=Deta(:,:,igrid)'*(Rhat(:,igrid).*(time>=times(sgrid)))/n;
        qgamhat(:,sgrid,igrid)=Dgam(:,:,igrid)'*(Rhat(:,igrid).*(time>=times(sgrid)))/n;
        qkq(:,sgrid)=qkq(:,sgrid)+(qgamhat(:,sgrid,igrid)-K(:,:,igrid)*qetahat(:,sgrid,igrid))*dt;
    end
end

for i=1:n
    Ti=time(i);
    Ci=cause(i);
    Mic=(Ti<=times)*(Ci==0)...
        -((Ti>=times).*dlogGc)*triu(ones(ngrid));
    dMic=Mic-[0,Mic(1:(ngrid-1))];
    y=sum(repmat(time,1,ngrid)>=repmat(times,n,1),1)/n;
    dMicyinv=dMic./y;
    
    for igrid=1:ngrid
        Phietaihat(:,i,igrid)=qetahat(:,:,igrid)*dMicyinv';
    end
    PhiGihat(:,i)=qkq*dMicyinv';
end

Igaminvhat=n*Igaminv;
Betaihat=zeros(px,n,ngrid);
Bgamihat=zeros(pz,n,ngrid);
BKBihat=zeros(pz,n,ngrid);

SWgamihat2=zeros(pz,pz);
for i=1:n
    for igrid=1:ngrid
        Betaihat(:,i,igrid)=Deta(i,:,igrid)'*(Rhat(i,igrid)-F1hat(i,igrid));
        Bgamihat(:,i,igrid)=Dgam(i,:,igrid)'*(Rhat(i,igrid)-F1hat(i,igrid));
        
        BKBihat(:,i,igrid)=Bgamihat(:,i,igrid)-K(:,:,igrid)*Betaihat(:,i,igrid);
    end
    BKBihatint=sum(BKBihat,3)*dt;
    Wgamihat=BKBihatint+PhiGihat;
    SWgamihat2=SWgamihat2+Wgamihat(:,i)*Wgamihat(:,i)';
end

Sigmagam=Igaminvhat*SWgamihat2*Igaminvhat/n/n;
sigmagam=sqrt(diag(Sigmagam));

Wetaihat=zeros(px,n,ngrid);
SWetaihat2=zeros(px,px,ngrid);
Sigmaeta=zeros(px,px,ngrid);
sigmaeta=zeros(px,ngrid);
for igrid=1:ngrid
    for i=1:n
        Wetaihat(:,i,igrid)=Betaihat(:,i,igrid)...
            +Phietaihat(:,i,igrid)...
            -Cgameta(:,:,igrid)'/n*Igaminvhat*Wgamihat(:,i);
        SWetaihat2(:,:,igrid)=SWetaihat2(:,:,igrid)...
            +Wetaihat(:,i,igrid)*Wetaihat(:,i,igrid)';
    end
    Sigmaeta(:,:,igrid)=pinv(Cetaeta(:,:,igrid))...
        *SWetaihat2(:,:,igrid)*pinv(Cetaeta(:,:,igrid));
    sigmaeta(:,igrid)=sqrt(diag(Sigmaeta(:,:,igrid)));
end







