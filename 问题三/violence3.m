clc;
clear all;

%%问题一求解
max_xy=30;   %3600/120=30
totalkm=40.168;   %大交路公里数
OD=load("OD.mat");
OD=table2array(OD.OD);
DM=load("DM.mat");
DM=table2array(DM.Duanmian);
dkm=load("km.mat");
dkm=table2array(dkm.Untitled);
runT=load("runT.mat");
runT=table2array(runT.runT);

%%评价变量矩阵，下标为可行方案的编号
proT=zeros(1,200);   %总等待时间
pros=zeros(1,200);   %总列车数
prokm=zeros(1,200);   %总公里数
%%自变量矩阵
prox=zeros(1,200);   %大交路列车数
proE=zeros(1,200);   %小交路终点站
proB=zeros(1,200);   %小交路起点站 
proStop=zeros(1,29,200);

proi=0;   %可行方案计数，兼最优解标记

%%评价变量最值初始化
max_T=0;
max_s=0;
max_km=0;
min_T=100000000;
min_s=100;
min_km=100000000;

cstop=[1,2,5,8,10,14,17,18,21,22,25,26,27,30];   %可以作为起始站终点站的车站

for Bi=1:length(cstop)-1   %小交路起始站begin
    for Ei=Bi+1:length(cstop)   %小交路终点站end
        if cstop(1,Ei)-cstop(1,Bi)<2   %小交路区间通过车站数量约束
            continue;
        end
        if cstop(1,Ei)-cstop(1,Bi)>23
            continue;
        end
        B=cstop(1,Bi);
        E=cstop(1,Ei);

        for x=1:max_xy   %大交路
            for y=1:max_xy-x   %小交路
                flag=1;   %标志，flag==0时受到约束跳出循环
                if y<10-x
                    continue;
                end
                xdown=zeros(1,30);   %只能坐大交路来到此站下车的人数
                xup=zeros(1,30);   %只能坐大交路在此站上车的人数
                xydown=zeros(1,30);   %可以坐大交路或者小交路来到此站下车的人数
                xyup=zeros(1,30);   %可以坐大交路或者小交路在此站上车的人数
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
                interT=ceil(3600/(x+y));   %发车间隔
                stop=zeros(1,29);          %记录可以不停的车站
                fx=zeros(1,29);   %断面列车人数约束，需要每站均小于1860
                fy=zeros(1,29);
                litkm=0;   %小交路公里数
                DMy=zeros(29,1);
                for i=1:29
                    if i>=B && i<E
                        litkm=litkm+dkm(i,1);
                        for j=B:i
                            for k=i+1:E
                                DMy(i,1)=DMy(i,1)+OD(j,k);
                                fy(1,i)=max(0,DMy(i,1)-1860*y);
                            end
                        end
                        fx(1,i)=ceil((DM(i,1)-DMy(i,1)+fy(1,i))/x);   %能坐小交路直达的都坐小交路，坐不下的再坐大交路
                    else 
                        fx(1,i)=ceil(DM(i,1)/x);
                    end
                    if fx(1,i)>1860   %断面列车人数约束
                        flag=0;
                        break;
                    end
                end
                if flag == 0
                    continue;
                end
                %%
                stop=zeros(1,29); 
                for i=B+1:E-1
                    stop(1,i)=(x+y)-ceil(DM(i,1)/1860);
                end
                stopmin=100;
                for i=1:29
                    if stop(1,i)<stopmin&&stop(1,i)>1
                        stopmin=stop(1,i);
                    end
                end
                if stopmin<100
                   alpha=stopmin;
                else
                    alpha=0;
                end
                deltaT=zeros(1,30);
                for i=2:29
                    if stop(1,i)>1
                        deltaT(1,i)=(xup(1,i)/x+xdown(1,i)/x+xydown(1,i)/(x+y-alpha)+xyup(1,i)/(x+y-alpha))*0.04;
                    else 
                        deltaT(1,i)=(xup(1,i)/x+xdown(1,i)/x+xydown(1,i)/(x+y)+xyup(1,i)/(x+y))*0.04;
                    end
                end
                %%
                for i=1:30
                    if deltaT(1,i)<20   %停站时间约束
                        deltaT(1,i)=20;
                    else
                        if deltaT(1,i)>120
                            flag=0;
                            break;
                        else
                            if i>1
                                if runT(i-1,1)>interT&&interT-deltaT(1,i)<108   %追踪时间约束
                                    flag=0;
                                    break;
                                end
                            end
                        end
                    end
                end
                if flag == 0
                    continue;
                end
                deltaT(1,1)=0;   %起点站停站时间不计入等待时间

                X=x;
                Y=y;
                if x<y   %使X>Y恒成立
                    X=y;
                    Y=x;
                end
                Q=floor(X/Y);
                remxy=rem(X,Y);
                if remxy>Q   %交替开行约束
                    flag=0;
                end
                if flag == 0
                    continue;
                end
                
                proi=proi+1;
                prokm(1,proi)=x*totalkm+y*litkm;   %总公里数
                totalT=0;   %总等待时间
                for i=1:29
                    for j=i+1:30
                        for di=i:j-1
                            totalT=totalT+OD(i,j)*deltaT(1,di);
                        end
                        if i>=B && j<=E
                            totalT=totalT+OD(i,j)*interT/2*((x+y)/(x+y-alpha));
                        else
                            totalT=totalT+OD(i,j)*interT*((x+y)/(x+y-alpha))*(x+y)/x/2;
                        end
                    end
                end
                
                %%数据存储
                prox(1,proi)=x;
                proE(1,proi)=E;
                proB(1,proi)=B;
                proA(1,proi)=alpha;
                pros(1,proi)=x+y;
                proT(1,proi)=totalT;
                proStop(:,:,proi)=stop;
                %%求最值，以进行归一化和求最优解
                if pros(1,proi)<min_s
                    min_s=pros(1,proi);
                end
                if pros(1,proi)>max_s
                    max_s=pros(1,proi);
                end
                if prokm(1,proi)<min_km
                    min_km=prokm(1,proi);
                end
                if prokm(1,proi)>max_km
                    max_km=prokm(1,proi);
                end
                if proT(1,proi)<min_T
                    min_T=proT(1,proi);
                end
                if proT(1,proi)>max_T
                    max_T=proT(1,proi);
                end
             end
        end
    end
end
%%求变化范围
avers=(max_s-min_s);
averT=(max_T-min_T);
averkm=(max_km-min_km);
%归一化平衡量纲

s_one=(pros-min_s)/avers;
T_one=(proT-min_T)/averT;
km_one=(prokm-min_km)/averkm;

profit = (pros-min_s)/avers+2*((proT-min_T)/averT)+(prokm-min_km)/averkm;
minpro=3;   %最优解值
xi=0;    %最优解标记
for i=1:length(profit)
    if profit(1,i)<minpro
        minpro=profit(1,i);
        xi=i;
    end
end
proy=pros-prox;
%导出最优解
minx=prox(1,xi);
miny=proy(1,xi);
minE=proE(1,xi);
minB=proB(1,xi);
mina=proA(1,xi);
minstop=proStop(:,:,xi);
proii=linspace(1,length(profit),length(profit));   %下标坐标


minpro
s_one(1,xi)
T_one(1,xi)
km_one(1,xi)
