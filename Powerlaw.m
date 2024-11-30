clc;close all
clear all
H=0.2;
M=22;
N=6;
x0=0;
y0=0;
X=zeros(N,M);
Y=zeros(N,M);
deltax=1/(M-2);
deltay=0.2/(N-2);
x0=deltax/2;
y0=deltay/2;
for i=1:N
    for j=2:M-1
        X(i,j)=x0+(j-2)*deltax;
    end
 
end
   X(:,M)=X(:,M-1)+deltax/2;
   for j=1:M
    for i=2:N-1
        Y(i,j)=y0+(i-2)*deltay;
    end
    
end
   Y(N,:)=Y(N-1,:)+deltay/2;
   for i=1:N
     j=1:M
     U(i,j)=1.5-1.5*((2*Y(i,j)/H)-1).^2
   end
kcp=0.02;
D=kcp;
ro=1;
omega=0.5;
 T=zeros(N,M);
 T(:,:)=10 ;
 T(1,:)=100;
 T(N,:)=100;
 T(:,1)=0;
 Tnew=T;
 error=1000;
 eps=.001;
  T(:,M)=T(:,M-1)
  
         
         while(error>=eps)
 
 for i=2:N-1
     for j=2:M-1
         deltaxe=X(i,j+1)-X(i,j);
         deltaxw=X(i,j)-X(i,j-1);
         deltayn=Y(i+1,j)-Y(i,j);
         deltays=Y(i,j)-Y(i-1,j);
         Ae=Y(i+1,j+1)-Y(i,j+1);
         Aw=Y(i+1,j)-Y(i,j);
         An=X(i+1,j+1)-X(i+1,j);
         As=X(i,j+1)-X(i,j);
         Fe=ro*U(i,j)*Ae;
         Fw=ro*U(i,j)*Aw;
         Fn=0;
         Fs=0;
         sp=0.;
         Su1=0; Su2=0;Su3=0;Su4=0;
         aS=D;
         aN=D;
        pee=Fe/D;
        pew=Fw/D;
        betae=((1-0.1*pee)^5);
        betaw=((1-0.1*pew)^5);
         aE=(D*betae)-Fe+D;
         aW=(D*betaw)+Fw+D;
         
        if j==2
        sp=-aW;     
        Su3=aW*T(i,j-1);
        aW=0;
        
          if (i==N-1)
             aN=2*aN;
             Su2=aN*T(i+1,j);
             sp=-aN;
             aN=0;
        end
        if (i==2)
            aS=2*aS;
             Su1=aS*T(i-1,j);
             sp=-aS;
             aS=0;
        end
        end
                 
        b=Su1+Su2+Su3+Su4;

         aP=aE+aW+aN+aS-sp+Fe-Fw;
         Tnew(i,j)=omega*((aE*T(i,j+1)+aW*T(i,j-1)+aN*T(i+1,j)+aS*T(i-1,j)+b)/aP)+(1-omega)*T(i,j) ;       
        end
        
        end
   Tnew(:,M)=Tnew(:,M-1);
error=norm((Tnew-T),2)
T=Tnew;
 end

 
contourf(X, Y, T, 20, 'LineColor', 'none');
colorbar;
title('Temperature Distribution in the Channel Power-law Method');
xlabel('x');
ylabel('y');
