% this program solve 2d diffusion equation
clc;close all
clear all
gama1=2000;%=T
M=32;
N=22;
x0=0;
y0=0;
X=zeros(N,M);
Y=zeros(N,M);
deltax=600/(M-2);
deltay=400/(N-2);
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

 B=50;
 h=zeros(N,M);
 S=zeros(N,M); % source term
 h(:,:)=0;
 h(N,:)=100;
 hnew=h;
 gama=zeros(N,M);
 gama(:,:)=gama1;
 omega=1;
 error=1000;
 eps=.001;
 qw=0.1;
 qe=0.0;
 qs=0.0;
 
 while(error>=eps)
 
 for j=2:M-1
     for i=2:N-1
         gamae=(gama(i,j+1)+gama(i,j))/2;
         gamaw=(gama(i,j-1)+gama(i,j))/2;
         gamas=(gama(i-1,j)+gama(i,j))/2;
         gaman=(gama(i+1,j)+gama(i,j))/2;
         deltaxe=X(i,j+1)-X(i,j);
         deltaxw=X(i,j)-X(i,j-1);
         deltayn=Y(i+1,j)-Y(i,j);
         deltays=Y(i,j)-Y(i-1,j);
         Ae=(Y(i+1,j+1)-Y(i,j+1))*B;
         Aw=(Y(i+1,j)-Y(i,j))*B;
         An=(X(i+1,j+1)-X(i+1,j))*B;
         As=(X(i,j+1)-X(i,j))*B;
                 
         aE=gamae*Ae/deltaxe;
         aW=gamaw*Aw/deltaxw;
         aN=gaman*An/deltayn;
         aS=gamas*As/deltays;
         SP=0.;
         b=0.;
         Su=S(i,j)*deltax*deltay;Su1=0; Su2=0;Su3=0;Su4=0;SR=0;
        if (i==N-1)
             aN=aN;
             Su2=aN*h(i+1,j);
             SP=-aN;
            aN=0;
        end
        if(j==2)
            aW=0.;
            Su3=qw*B*deltay;
        end
        if(j==M-1)
            aE=0.;
            Su4=qe*B*deltay;
        end
        if (i==2)
             aS=0.;
             Su1=qs*B*deltax;
        end

         aP=aE+aW+aN+aS-SP;
           
             
         hnew(i,j)=omega*(aE*h(i,j+1)+aW*h(i,j-1)+aS*h(i-1,j)+aN*h(i+1,j)+Su1+Su2+Su3+Su4+Su+SR)/aP+(1-omega)*h(i,j);
     end
  
 end
hnew(1,:)=hnew(2,:);hnew(:,M)=hnew(:,M-1);hnew(:,1)=hnew(:,2)+qw*2*deltax/(2000);hnew(N,:)=100;
error=norm((hnew-h),2)
h=hnew;
 end
 contourf(X, Y, h, 20, 'LineColor', 'none');
colorbar;
title('water head  Distribution  Uniform Grid');
xlabel('x');
ylabel('y');


