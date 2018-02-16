clear;
clc;
L=60;
W=40;
Vo=5;
% (a)

k=0:0.5:L;
V=@(k)Vo-Vo/L*k;

plot(k,V(k));
title('part(a)');
xlabel('x');
ylabel('V');


% discrete=@(n) 1/n*cosh(n*pi*x/L)/cosh(n*pi*W/L)*sin(n*pi*y/L);
% V=@(x,y) 4*Vo/pi*symsum(discrete(2*k+1),k,1,Inf);

% (b)
G=zeros(L,W);
[X,Y]=meshgrid(1:1:W,1:1:L);

for x=1:L
    for y=1:W
%         if x==1 || x==L
%             voltage(x,y)=Vo;
%         elseif y==1 || y==W
%             voltage(x,y)=0;
%         else
        G(x,y)=0;
        for n=1:2:191
            G(x,y)=G(x,y)+4*Vo/pi*1/n*cosh(n*pi*x/L)/cosh(n*pi*W/L)*sin(n*pi*y/L);
        end
%         end
    end
surf(X,Y,G);
title('part(b)');
xlabel('X');
ylabel('Y');

pause(0.1);
end

