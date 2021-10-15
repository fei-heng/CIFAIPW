function  [etahatall,sigmaeta,gamhat,sigmagam,Fhat,sigmaF]=esta(n,px,pz,...
    etahatall0,gamhat0,times,ngrid,dt,Gctimes,...
    time,X,Z,A,Auxind,Gcx,cause,w_est,...
    rselec,shat,dpi,Itheta,zetai,ptheta,...
    nF,pd)

delta=(cause==1);
Z(isnan(Z))=0;
z1=Z(:,1);
z2=Z(:,2);

ind_case=(cause==1);
ind_ntype=(cause==2);
ind_censor=(cause==0);

ind_case_ava=(cause==1)&(rselec==1);
ind_ntype_ava=(cause==2)&(rselec==1);
ind_censor_ava=(cause==0)&(rselec==1);

if Auxind==0
    J12=[X,z1,log(time)];
    J0=[X,z1];
else
    J12=[X,z1,A,log(time)];
    J0=[X,z1,A];
end
J12_case=J12(ind_case,:);
n1=sum(ind_case);
J12_ntype=J12(ind_ntype,:);
n2=sum(ind_ntype);
J0_censor=J0(ind_censor,:);
n3=sum(ind_censor);
J12_case_ava=J12(ind_case_ava,:);
n1_ava=sum(ind_case_ava);
J12_ntype_ava=J12(ind_ntype_ava,:);
n2_ava=sum(ind_ntype_ava);
J0_censor_ava=J0(ind_censor_ava,:);
n3_ava=sum(ind_censor_ava);

pnu1=size(J12,2);
pnu2=size(J12,2);
pnu3=size(J0,2);


F1hat=zeros(n,ngrid);
Rhat=zeros(n,ngrid);
Deta=zeros(n,px,ngrid);
Dgam=zeros(n,pz,ngrid);
%Wt=zeros(n,n,ngrid);
%Psihat=zeros(n,n,ngrid);
%Wt=eye(n);
Psihat=diag(w_est);
OnePsihat=diag(1-w_est);

Beta=zeros(px,ngrid);
Bgam=zeros(pz,ngrid);
BKB=zeros(pz,ngrid);

Cetaeta=zeros(px,px,ngrid);
Cgameta=zeros(pz,px,ngrid);
Cgamgam=zeros(pz,pz,ngrid);

K=zeros(pz,px,ngrid);
CKC=zeros(pz,pz,ngrid);

% esta function
etahatall=etahatall0;
gamhat=gamhat0;
error = 1;
iter = 1;
while error>0.00001 && iter<30
    %iter
    eeta=zeros(n,px,ngrid);
    egam=zeros(n,pz,ngrid);
    Cetaeta_a=zeros(px,px,ngrid);
    Cgameta_a=zeros(pz,px,ngrid);
    Cgamgam_a=zeros(pz,pz,ngrid);
    
    EDeta=zeros(n,px,ngrid);
    EDgam=zeros(n,pz,ngrid);
    
    eg1z1=exp(z1*gamhat(1));
    eg2z2=exp(z2*gamhat(2));
    
    nu2a=zeros(pnu2,ngrid);
    nu2b=zeros(pnu2,ngrid);
    nu2c=zeros(pnu2,ngrid);
    nu2d=zeros(pnu2,ngrid);
    
    nu3a=zeros(pnu3,ngrid);
    nu3b=zeros(pnu3,ngrid);
    nu3c=zeros(pnu3,ngrid);
    nu3d=zeros(pnu3,ngrid);
    
    %         a_ntype_ava=zeros(n2_ava,ngrid);
    %         b_ntype_ava=zeros(n2_ava,ngrid);
    %         c_ntype_ava=zeros(n2_ava,ngrid);
    %         d_ntype_ava=zeros(n2_ava,ngrid);
    %
    %         a_censor_ava=zeros(n3_ava,ngrid);
    %         b_censor_ava=zeros(n3_ava,ngrid);
    %         c_censor_ava=zeros(n3_ava,ngrid);
    %         d_censor_ava=zeros(n3_ava,ngrid);
    
    a=zeros(n,ngrid);
    b=zeros(n,ngrid);
    c=zeros(n,ngrid);
    d=zeros(n,ngrid);
    
    for igrid=1:ngrid
        %igrid
        t=times(igrid);
        etahat=etahatall(igrid,:);
        Xeta=X*etahat';
        
        Rhat(:,igrid)=(time<=t).*delta./Gcx;
        F1hat(:,igrid)=1-exp(-Xeta.*exp(Z*gamhat'));
        Deta(:,:,igrid)=repmat(exp(-Xeta.*exp(Z*gamhat')).*exp(Z*gamhat'),1,px).*X;
        Dgam(:,:,igrid)=repmat(exp(-Xeta.*exp(Z*gamhat')).*Xeta.*exp(Z*gamhat'),1,pz).*Z;
        %Wt(:,:,igrid)=eye(n);
        %Psihat(:,:,igrid)=diag(w_est);
        
        %C1=X.*repmat(eg1z1.*(Rhat(:,igrid)-1),1,px);
        %C2=exp(-(X*etahat').*eg1z1);
        V=exp(-Xeta.*eg1z1).^eg2z2;
        Ea=V.*eg2z2;
        Eb=Ea.*z2;
        Ec=Ea.*V;
        Ed=Ec.*z2;
        Ee=Ed.*eg2z2;
        Ef=Ec.*eg2z2;
        Eg=Ee.*z2;
        
        if sum(isnan(Eg)+isinf(Eg)+isnan(eg1z1)+isinf(eg1z1))>0
            etahatall=NaN;
            gamhat=NaN;
            sigmaeta=NaN;
            sigmagam=NaN;
            Fhat=NaN;
            sigmaF=NaN;
            return
        end
        
        a(:,igrid)=Ea;
        b(:,igrid)=Eb;
        c(:,igrid)=Ec;
        d(:,igrid)=Ed;
        
        
        % add nu1a-g later
        %a_ntype_ava(:,igrid)=Ea(ind_ntype_ava);
        nu2a(:,igrid)=regress(Ea(ind_ntype_ava),J12_ntype_ava);
        Ea(ind_ntype)=J12_ntype*nu2a(:,igrid);
        %b_ntype_ava(:,igrid)=Eb(ind_ntype_ava);
        nu2b(:,igrid)=regress(Eb(ind_ntype_ava),J12_ntype_ava);
        Eb(ind_ntype)=J12_ntype*nu2b(:,igrid);
        %c_ntype_ava(:,igrid)=Ec(ind_ntype_ava);
        nu2c(:,igrid)=regress(Ec(ind_ntype_ava),J12_ntype_ava);
        Ec(ind_ntype)=J12_ntype*nu2c(:,igrid);
        %d_ntype_ava(:,igrid)=Ed(ind_ntype_ava);
        nu2d(:,igrid)=regress(Ed(ind_ntype_ava),J12_ntype_ava);
        Ed(ind_ntype)=J12_ntype*nu2d(:,igrid);
        nu2e=regress(Ee(ind_ntype_ava),J12_ntype_ava);
        Ee(ind_ntype)=J12_ntype*nu2e;
        nu2f=regress(Ef(ind_ntype_ava),J12_ntype_ava);
        Ef(ind_ntype)=J12_ntype*nu2f;
        nu2g=regress(Eg(ind_ntype_ava),J12_ntype_ava);
        Eg(ind_ntype)=J12_ntype*nu2g;
        
        %a_censor_ava(:,igrid)=Ea(ind_censor_ava);
        nu3a(:,igrid)=regress(Ea(ind_censor_ava),J0_censor_ava);
        Ea(ind_censor)=J0_censor*nu3a(:,igrid);
        %b_censor_ava(:,igrid)=Eb(ind_censor_ava);
        nu3b(:,igrid)=regress(Eb(ind_censor_ava),J0_censor_ava);
        Eb(ind_censor)=J0_censor*nu3b(:,igrid);
        %c_censor_ava(:,igrid)=Ec(ind_censor_ava);
        nu3c(:,igrid)=regress(Ec(ind_censor_ava),J0_censor_ava);
        Ec(ind_censor)=J0_censor*nu3c(:,igrid);
        %d_censor_ava(:,igrid)=Ed(ind_censor_ava);
        nu3d(:,igrid)=regress(Ed(ind_censor_ava),J0_censor_ava);
        Ed(ind_censor)=J0_censor*nu3d(:,igrid);
        nu3e=regress(Ee(ind_censor_ava),J0_censor_ava);
        Ee(ind_censor)=J0_censor*nu3e;
        nu3f=regress(Ef(ind_censor_ava),J0_censor_ava);
        Ef(ind_censor)=J0_censor*nu3f;
        nu3g=regress(Eg(ind_censor_ava),J0_censor_ava);
        Eg(ind_censor)=J0_censor*nu3g;
       
        
        eeta(:,:,igrid)=X.*repmat(eg1z1.*(Ea.*(Rhat(:,igrid)-1)+Ec),1,px);
        egam(:,:,igrid)=[z1.*(Ea.*(Rhat(:,igrid)-1)+Ec),...
            (Eb.*(Rhat(:,igrid)-1)+Ed)].*repmat(Xeta.*eg1z1,1,pz);
        
        EDeta(:,:,igrid)=X.*repmat(Ea.*eg1z1,1,px);
        EDgam(:,:,igrid)=[z1.*Ea,Eb].*repmat(Xeta.*eg1z1,1,pz);
        
        for i=1:n
            Cetaeta_a(:,:,igrid)=Cetaeta_a(:,:,igrid)...
                +X(i,:)'*X(i,:)*(1-w_est(i))*Ef(i)*eg1z1(i)^2;
            Cgameta_a(:,:,igrid)=Cgameta_a(:,:,igrid)...
                +[Ef(i)*z1(i);Ee(i)]*X(i,:)*(1-w_est(i))*Xeta(i)*eg1z1(i)^2;
            Cgamgam_a(:,:,igrid)=Cgamgam_a(:,:,igrid)...
                +[Ef(i)*z1(i)^2,Ee(i)*z1(i);Ee(i)*z1(i),Eg(i)]*(1-w_est(i))*(Xeta(i)*eg1z1(i))^2;
        end
        
        Beta(:,igrid)=Deta(:,:,igrid)'*Psihat*(Rhat(:,igrid)-F1hat(:,igrid))+eeta(:,:,igrid)'*(1-w_est);
        Bgam(:,igrid)=Dgam(:,:,igrid)'*Psihat*(Rhat(:,igrid)-F1hat(:,igrid))+egam(:,:,igrid)'*(1-w_est);
        
        
        Cetaeta(:,:,igrid)=Deta(:,:,igrid)'*Psihat*Deta(:,:,igrid)+Cetaeta_a(:,:,igrid);
        Cgameta(:,:,igrid)=Dgam(:,:,igrid)'*Psihat*Deta(:,:,igrid)+Cgameta_a(:,:,igrid);
        Cgamgam(:,:,igrid)=Dgam(:,:,igrid)'*Psihat*Dgam(:,:,igrid)+Cgamgam_a(:,:,igrid);
        
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
    etahatall=NaN;
    gamhat=NaN;
    sigmaeta=NaN;
    sigmagam=NaN;
    Fhat=NaN;
    sigmaF=NaN;
end

if any(isnan(etahatall))
    return
end

% dmuetadnu1
% dmugamdnu2
gnu1t_nu2a=zeros(px,pnu2,ngrid);
gnu1t_nu3a=zeros(px,pnu3,ngrid);
gnu1t_nu2c=zeros(px,pnu2,ngrid);
gnu1t_nu3c=zeros(px,pnu3,ngrid);
gnu1t=zeros(px,(pnu2+pnu3)*2,ngrid);

gnu2t_nu2a=zeros(pnu2,ngrid);
gnu2t_nu3a=zeros(pnu3,ngrid);
gnu2t_nu2b=zeros(pnu2,ngrid);
gnu2t_nu3b=zeros(pnu3,ngrid);
gnu2t_nu2c=zeros(pnu2,ngrid);
gnu2t_nu3c=zeros(pnu3,ngrid);
gnu2t_nu2d=zeros(pnu2,ngrid);
gnu2t_nu3d=zeros(pnu3,ngrid);
gnu2t=zeros(pz,(pnu2+pnu3)*4,ngrid);

for igrid=1:ngrid
    t=times(igrid);
    etahat=etahatall(igrid,:);
    Xeta=X*etahat';
    
    for i=1:n
        if ind_ntype(i)
            gnu1t_nu2a(:,:,igrid)=gnu1t_nu2a(:,:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*X(i,:)'*J12(i,:);
            gnu1t_nu2c(:,:,igrid)=gnu1t_nu2c(:,:,igrid)...
                +(1-w_est(i))*eg1z1(i)*X(i,:)'*J12(i,:);
            gnu2t_nu2a(:,igrid)=gnu2t_nu2a(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*Xeta(i)*z1(i)*J12(i,:)';
            gnu2t_nu2b(:,igrid)=gnu2t_nu2b(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*Xeta(i)*J12(i,:)';
            gnu2t_nu2c(:,igrid)=gnu2t_nu2c(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*Xeta(i)*z1(i)*J12(i,:)';
            gnu2t_nu2d(:,igrid)=gnu2t_nu2d(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*Xeta(i)*J12(i,:)';
        elseif ind_censor(i)
            gnu1t_nu3a(:,:,igrid)=gnu1t_nu3a(:,:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*X(i,:)'*J0(i,:);
            gnu1t_nu3c(:,:,igrid)=gnu1t_nu3c(:,:,igrid)...
                +(1-w_est(i))*eg1z1(i)*X(i,:)'*J0(i,:);
            gnu2t_nu3a(:,igrid)=gnu2t_nu3a(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*Xeta(i)*z1(i)*J0(i,:)';
            gnu2t_nu3b(:,igrid)=gnu2t_nu3b(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*(Rhat(i)-1)*Xeta(i)*J0(i,:)';
            gnu2t_nu3c(:,igrid)=gnu2t_nu3c(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*Xeta(i)*z1(i)*J0(i,:)';
            gnu2t_nu3d(:,igrid)=gnu2t_nu3d(:,igrid)...
                +(1-w_est(i))*eg1z1(i)*Xeta(i)*J0(i,:)';
        end
    end
    gnu1t(:,:,igrid)=[gnu1t_nu2a(:,:,igrid),gnu1t_nu3a(:,:,igrid),...
        gnu1t_nu2c(:,:,igrid),gnu1t_nu3c(:,:,igrid)]/n;
    gnu2t(:,:,igrid)=[gnu2t_nu2a(:,igrid)',gnu2t_nu3a(:,igrid)',zeros(1,pnu2+pnu3),...
        gnu2t_nu2c(:,igrid)',gnu2t_nu3c(:,igrid)',zeros(1,pnu2+pnu3);...
        zeros(1,pnu2+pnu3),gnu2t_nu2b(:,igrid)',gnu2t_nu3b(:,igrid)',...
        zeros(1,pnu2+pnu3),gnu2t_nu2d(:,igrid)',gnu2t_nu3d(:,igrid)']/n;
end

% I2 for nu2a-d
% I3 for nu3a-d
I2=zeros(pnu2,pnu2);
I3=zeros(pnu3,pnu3);
for i=1:n2_ava
    I2=I2+J12_ntype_ava(i,:)'*J12_ntype_ava(i,:)/n2_ava;
end
I2inv=pinv(I2);
for i=1:n3_ava
    I3=I3+J0_censor_ava(i,:)'*J0_censor_ava(i,:)/n3_ava;
end
I3inv=pinv(I3);
Inu1inv=blkdiag(I2inv,I3inv,I2inv,I3inv);
Inu2inv=blkdiag(I2inv,I3inv,I2inv,I3inv,I2inv,I3inv,I2inv,I3inv);

phinua2=zeros(pnu2,n,ngrid);
phinua3=zeros(pnu3,n,ngrid);
phinub2=zeros(pnu2,n,ngrid);
phinub3=zeros(pnu3,n,ngrid);
phinuc2=zeros(pnu2,n,ngrid);
phinuc3=zeros(pnu3,n,ngrid);
phinud2=zeros(pnu2,n,ngrid);
phinud3=zeros(pnu3,n,ngrid);
phinu1=zeros((pnu2+pnu3)*2,n,ngrid);
phinu2=zeros((pnu2+pnu3)*4,n,ngrid);
for igrid=1:ngrid
    for i=1:n
        if ind_ntype_ava(i)==1
            phinua2(:,i,igrid)=(a(i,igrid)-J12(i,:)*nu2a(:,igrid))*J12(i,:)'*n/n2_ava;
            phinub2(:,i,igrid)=(b(i,igrid)-J12(i,:)*nu2b(:,igrid))*J12(i,:)'*n/n2_ava;
            phinuc2(:,i,igrid)=(c(i,igrid)-J12(i,:)*nu2c(:,igrid))*J12(i,:)'*n/n2_ava;
            phinud2(:,i,igrid)=(d(i,igrid)-J12(i,:)*nu2d(:,igrid))*J12(i,:)'*n/n2_ava;
        elseif ind_censor_ava(i)==1
            phinua3(:,i,igrid)=(a(i,igrid)-J0(i,:)*nu3a(:,igrid))*J0(i,:)'*n/n3_ava;
            phinub3(:,i,igrid)=(b(i,igrid)-J0(i,:)*nu3b(:,igrid))*J0(i,:)'*n/n3_ava;
            phinuc3(:,i,igrid)=(c(i,igrid)-J0(i,:)*nu3c(:,igrid))*J0(i,:)'*n/n3_ava;
            phinud3(:,i,igrid)=(d(i,igrid)-J0(i,:)*nu3d(:,igrid))*J0(i,:)'*n/n3_ava;
        end
    end
    phinu1(:,:,igrid)=vertcat(phinua2(:,:,igrid),phinua3(:,:,igrid),...
        phinuc2(:,:,igrid),phinuc3(:,:,igrid));
    phinu2(:,:,igrid)=vertcat(phinua2(:,:,igrid),phinua3(:,:,igrid),...
        phinub2(:,:,igrid),phinub3(:,:,igrid),...
        phinuc2(:,:,igrid),phinuc3(:,:,igrid),...
        phinud2(:,:,igrid),phinud3(:,:,igrid));
end

% Phieta, PhiG
% geta, ggam, Ithetainv, zetai, Phipi
Phietaihat=zeros(px,n,ngrid);
PhiGihat=zeros(pz,n);
Phimu=zeros(pz,n);

qetahat=zeros(px,ngrid,ngrid);
qgamhat=zeros(pz,ngrid,ngrid);
qkq=zeros(pz,ngrid);

logGc=-log(Gctimes);
dlogGc=logGc-[0,logGc(1:(ngrid-1))];

for sgrid=1:ngrid
    for igrid=1:ngrid
        qetahat(:,sgrid,igrid)=(Deta(:,:,igrid)'*Psihat+EDeta(:,:,igrid)'*OnePsihat)*(Rhat(:,igrid).*(time>=times(sgrid)))/n;
        qgamhat(:,sgrid,igrid)=(Dgam(:,:,igrid)'*Psihat+EDgam(:,:,igrid)'*OnePsihat)*(Rhat(:,igrid).*(time>=times(sgrid)))/n;
        qkq(:,sgrid)=qkq(:,sgrid)+(qgamhat(:,sgrid,igrid)-K(:,:,igrid)*qetahat(:,sgrid,igrid))*dt;
    end
end

geta=zeros(px,ptheta,ngrid);
ggam=zeros(pz,ptheta);
Ithetainv=pinv(Itheta);
for i=1:n
    for igrid=1:ngrid
        geta(:,:,igrid)=geta(:,:,igrid)-rselec(i)*(Deta(i,:,igrid)'...
            *(Rhat(i,igrid)-F1hat(i,igrid))-eeta(i,:,igrid)')/shat(i)^2*dpi(:,i)'/n;
        ggam(:,:)=ggam(:,:)-rselec(i)*(Dgam(i,:,igrid)'...
            *(Rhat(i,igrid)-F1hat(i,igrid))-egam(i,:,igrid)')/shat(i)^2*dpi(:,i)'/n*dt;
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
        Phietaihat(:,i,igrid)=qetahat(:,:,igrid)*dMicyinv'...
            +geta(:,:,igrid)*Ithetainv*zetai(:,i)+gnu1t(:,:,igrid)*Inu1inv*phinu1(:,i,igrid);
        Phimu(:,i)=Phimu(:,i)+gnu2t(:,:,igrid)*Inu2inv*phinu2(:,i,igrid)*dt...
            -K(:,:,igrid)*gnu1t(:,:,igrid)*Inu1inv*phinu1(:,i,igrid)*dt;
    end
    PhiGihat(:,i)=qkq*dMicyinv';
end

kgeta=zeros(pz,ptheta);
for igrid=1:ngrid
    kgeta=kgeta+K(:,:,igrid)*geta(:,:,igrid)*dt;
end

Igaminvhat=n*Igaminv;
Betaihat=zeros(px,n,ngrid);
Bgamihat=zeros(pz,n,ngrid);
BKBihat=zeros(pz,n,ngrid);

SWgamihat2=zeros(pz,pz);
for i=1:n
    for igrid=1:ngrid
        Betaihat(:,i,igrid)=Deta(i,:,igrid)'*w_est(i)*(Rhat(i,igrid)-F1hat(i,igrid))+eeta(i,:,igrid)'*(1-w_est(i));
        Bgamihat(:,i,igrid)=Dgam(i,:,igrid)'*w_est(i)*(Rhat(i,igrid)-F1hat(i,igrid))+egam(i,:,igrid)'*(1-w_est(i));
        
        BKBihat(:,i,igrid)=Bgamihat(:,i,igrid)-K(:,:,igrid)*Betaihat(:,i,igrid);
    end
end
BKBihatint=sum(BKBihat,3)*dt;

Phipi=(ggam-kgeta)*Ithetainv*zetai;
Wgamihat=BKBihatint+PhiGihat+Phipi+Phimu;
for i=1:n
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

% %% cif
% Fhat = zeros(ngrid,nF);
% for iF=1:nF
%    xpd=pd(iF,1:px); 
%    zpd=pd(iF,px+1:end);
%    Fhat(:,iF)=cifest(etahatall,xpd,gamhat,zpd)';
% end
% 
% WFihat=zeros(n,ngrid,nF);
% sigmaF=zeros(ngrid,nF);
% for iF=1:nF
%     xpd=pd(iF,1:px);
%     zpd=pd(iF,px+1:end);
%     for igrid=1:ngrid
%         for i=1:n
%             WFihat(i,igrid,iF)=exp(zpd*gamhat')*xpd...
%                 *pinv(Cetaeta(:,:,igrid))*Wetaihat(:,i,igrid)*n...
%                 +etahatall(igrid,:)*xpd'*exp(zpd*gamhat')...
%                 *zpd*Igaminvhat*Wgamihat(:,i);
%         end
%         sigmaF(igrid,iF)=(1-Fhat(igrid,iF))*sqrt(sum(WFihat(:,igrid,iF).*WFihat(:,igrid,iF)))/n;
%     end
% end







