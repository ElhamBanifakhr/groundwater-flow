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
G=[0	50	70	85	102	117	139	164	187	287	302	319	343	371	386	407	417	428	438	449	461	474	488	503	519	536	554	573	580	591	595	600]
O=[0
52
56
62
70
80
92
106
122
140
160
182
206
232
260
290
322
329
355
362
387
400
]

for i=1:N
    for j=1:M
        
       X(:,j)=G(1,j)
       Y(i,:)=O(i,1)
    end
 

end
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
         Su1=0; Su2=0;Su3=0;Su4=0;SR=0;
        if (i==N-1)
             aN=aN;
             Su2=aN*h(i+1,j);
             SP=-aN;
            aN=0;
        end
        if(j==2)
            aW=0.;
            Su3=qw*B*(deltayn+deltays)/2;
        end
        if(j==M-1)
            aE=0.;
            Su4=0;
        end
        if (i==2)
             aS=0.;
             Su1=0;
        end

         aP=aE+aW+aN+aS-SP;
           
             
         hnew(i,j)=omega*(aE*h(i,j+1)+aW*h(i,j-1)+aS*h(i-1,j)+aN*h(i+1,j)+Su1+Su2+Su3+Su4)/aP+(1-omega)*h(i,j);
     end
  
 end
hnew(1,:)=hnew(2,:);hnew(:,M)=hnew(:,M-1);hnew(:,1)=hnew(:,2)+qw*2*50/2000;hnew(N,:)=100;
error=norm((hnew-h),2)
h=hnew;
 end
 contourf(X, Y, h, 20, 'LineColor', 'none');
colorbar;
title('water head  Distribution  Uniform Grid');
xlabel('x');
ylabel('y');
