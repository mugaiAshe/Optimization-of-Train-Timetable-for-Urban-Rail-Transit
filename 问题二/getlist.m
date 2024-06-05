clc;
clear all;

totalkm=40.168;   %大交路公里数
OD=load("OD.mat");
OD=table2array(OD.OD);
DM=load("DM.mat");
DM=table2array(DM.Duanmian);
dkm=load("km.mat");
dkm=table2array(dkm.Untitled);
runT=load("runT.mat");
runT=table2array(runT.runT);

B=10;
E=17;
x=19;   %大交路
y=6;   %小交路
xdown=linspace(0,0,30);   %只能坐大交路来到此站下车的人数
xup=linspace(0,0,30);   %只能坐大交路在此站上车的人数
xydown=linspace(0,0,30);   %可以坐大交路或者小交路来到此站下车的人数
xyup=linspace(0,0,30);   %可以坐大交路或者小交路在此站上车的人数
for i=1:30
    for j=i:30
        if i>=B&&j<=E
            xyup(1,i)=xyup(1,i)+OD(i,j);
        else
            xup(1,i)=OD(i,j)+xup(1,i);
        end
    end
    for j=1:i
        if j>=B&&i<=E
            xydown(1,i)=xydown(1,i)+OD(j,i);
        else
            xdown(1,i)=OD(j,i)+xdown(1,i);
        end
    end
end
deltaT=(xup/x+xdown/x+xydown/(x+y)+xyup/(x+y))*0.04;   %停站时间,先根据总人数计算总停车时间
interT=ceil(3600/(x+y));   %发车间隔
for i=1:30
    if deltaT(1,i)<20   %停站时间约束
        deltaT(1,i)=20;
    end
end
deltaT(1,1)=0;

nonT=zeros(58,x+y);   %初始化用时间值0矩阵
nowT=duration(nonT,nonT,nonT);   %输出的时刻表初始化
beginT=hours(7);   %起始时刻
endT=duration(8,0,0);
num=0;   %车次

remx=0;
remy=0;
if x<y   %使X>Y恒成立
    X=1;
    Y=floor(y/x);
    remy=rem(y,x);
else
    Y=1;
    X=floor(x/y);
    remx=rem(x,y);
end

if remy>0
    for p=1:remy
        num=num+1;
        flag=1;
        for k=2*B-1:2*(E-1)
            if flag == 0
                nowT(k,num)=duration(0,0,0);
                continue;
            end
            for r=k:2*(E-1)
                if k==2*B-1
                    nowT(r,num)=nowT(r,num)+beginT+seconds(interT*(num-1));
                else
                    if rem(k,2)==1
                        nowT(r,num)=nowT(r,num)+seconds(deltaT(1,(k+1)/2));
                    else
                        nowT(r,num)=nowT(r,num)+seconds(runT(k/2,1));
                    end
                end
                if nowT(r,num)>endT
                    nowT(r,num)=endT;
                    flag=0;
                    break;
                end
            end
        end
    end
end

for i=1:(x-remx)/X
    for t=1:X
        num=num+1;
        flag=1;
        for k=1:58
            if flag == 0
                nowT(k,num)=duration(0,0,0);
                continue;
            end
            for r=k:58
                if k==1
                    nowT(r,num)=nowT(r,num)+beginT+seconds(interT*(num-1));
                else
                    if rem(k,2)==1
                        nowT(r,num)=nowT(r,num)+seconds(deltaT(1,(k+1)/2));
                    else
                        nowT(r,num)=nowT(r,num)+seconds(runT(k/2,1));
                    end
                end
                if nowT(r,num)>endT
                    nowT(r,num)=endT;
                    flag=0;
                    break;
                end
            end
        end
    end
    for u=1:Y
        num=num+1;
        flag=1;
        for k=2*B-1:2*(E-1)
            if flag == 0
                nowT(k,num)=duration(0,0,0);
                continue;
            end
            for r=k:2*(E-1)
                if k==2*B-1
                    nowT(r,num)=nowT(r,num)+beginT+seconds(interT*(num-1));
                else
                    if rem(k,2)==1
                        nowT(r,num)=nowT(r,num)+seconds(deltaT(1,(k+1)/2));
                    else
                        nowT(r,num)=nowT(r,num)+seconds(runT(k/2,1));
                    end
                end
                if nowT(r,num)>endT
                    nowT(r,num)=endT;
                    flag=0;
                    break;
                end
            end
        end
    end
end

if remx>0
    for p=1:remx
        num=num+1;
        flag=1;
        for k=1:58
            if flag == 0
                nowT(k,num)=duration(0,0,0);
                continue;
            end
            for r=k:58
                if k==1
                    nowT(r,num)=nowT(r,num)+beginT+seconds(interT*(num-1));
                else
                    if rem(k,2)==1
                        nowT(r,num)=nowT(r,num)+seconds(deltaT(1,(k+1)/2));
                    else
                        nowT(r,num)=nowT(r,num)+seconds(runT(k/2,1));
                    end
                end
                if nowT(r,num)>endT
                    nowT(r,num)=endT;
                    flag=0;
                    break;
                end
            end
        end
    end
end

proix=linspace(1,58,58);   %下标坐标
proiy=linspace(1,25,25);   %下标坐标
grid on;
set(gca,'GridLineStyle',':');

endT=duration(8,0,0)
plot(nowT(proix,proiy),ceil((proix+1)/2))
xlim([beginT endT]);
