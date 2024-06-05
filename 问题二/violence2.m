clc;
clear all;

%%问题二求解
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
prof=zeros(1,200);   %总偏差变量
%%自变量矩阵
prox=zeros(1,200);   %大交路列车数
proE=zeros(1,200);   %小交路终点站
proB=zeros(1,200);   %小交路起点站

proi=0;   %可行方案计数，兼最优解标记

%%评价变量最值初始化
max_T=0;
max_s=0;
max_km=0;
max_f=0;
min_T=100000000;
min_s=100;
min_km=100000000;
min_f=100000000;

cstop=[1,2,5,8,10,14,17,18,21,22,25,26,27,30];   %可以作为起始站终点站的车站

for Bi=1:length(cstop)-1   %小交路起始站begin
    for Ei=Bi+1:length(cstop)   %小交路终点站end，同时进行区间通过车站数量约束
        if cstop(1,Ei)-cstop(1,Bi)<2
            continue;
        end
        if cstop(1,Ei)-cstop(1,Bi)>23
            continue;
        end
        B=cstop(1,Bi);
        E=cstop(1,Ei);

        for x=1:max_xy   %大交路数量x
            for y=1:max_xy-x   %小交路数量y
                flag=1;   %标志，flag==0时受到约束跳出循环
                if y<10-x
                    continue;
                end

                %%交替开行约束
                X=x;
                Y=y;
                if x<y   %使X>Y恒成立
                    X=y;
                    Y=x;
                end
                Q=floor(X/Y);
                remxy=rem(X,Y);
                if remxy>Q
                    flag=0;
                end
                if flag == 0
                    continue;
                end

                xdown=zeros(1,30);   %坐大交路来到此站下车的人数
                xup=zeros(1,30);   %坐大交路在此站上车的人数
                ydown=zeros(1,30);   %坐小交路来到此站下车的人数
                yup=zeros(1,30);   %坐小交路在此站上车的人数
                for i=1:30
                    for j=i:30
                        if i>=B&&j<=E
                            yup(1,i)=yup(1,i)+OD(i,j);
                        else
                            xup(1,i)=OD(i,j)+xup(1,i);
                        end
                    end
                    for j=1:i
                        if j>=B&&i<=E
                            ydown(1,i)=ydown(1,i)+OD(j,i);   %能坐小交路的全坐小交路
                        else
                            xdown(1,i)=OD(j,i)+xdown(1,i);
                        end
                    end
                end

                newOD=OD;   %refresh OD
                ODx_in_xy=zeros(30,30);   %newOD中被挤得只能乘坐小交路的OD客流
                max_px=x*1860;
                max_py=y*1860;
                personx=0;   %大交路车上现有乘客
                persony=0;   %小交路车上现有乘客
                deviaf=0;   %负偏差变量
                
                for i=1:29
                    if i>=B && i<=E  %在大小交路重合区间
                        %%下车
                        persony=persony-ydown(1,i);
                        personx=personx-xdown(1,i);
                        %%比较是否够载
                        if personx+persony+yup(1,i)+xup(1,i)<=max_px+max_py && personx+xup(1,i)<=max_px   %够载，没有偏差变量
                            if persony+yup(1,i)<=max_py
                                %%能坐小交路的全坐小交路
                                persony=persony+yup(1,i);
                                personx=personx+xup(1,i);
                            else
                                %%小交路满，大交路有空，部分可以坐小交路的人被挤去大交路
                                tempf=(yup(1,i)+persony)-max_py;
                                totempf=0;   %以totempf达到tempf为偏差已经完全分配到newOD中的标志
                                personx=personx+xup(1,i)+tempf;
                                persony=max_py;
                                for j=i+1:E
                                    pxy=floor(newOD(i,j)*tempf/yup(1,i));   %以tempf/yup(1,i)为比例向下取整减少客流量，并作相应记录
                                    totempf=totempf+pxy;
                                    xdown(1,j)=xdown(1,j)+pxy;
                                    ydown(1,j)=ydown(1,j)-pxy;
                                    ODx_in_xy(i,j)=pxy;
                                end
                                tempj=i;
                                while totempf<tempf && tempj<E   %将偏差的零头由近及远分配
                                    tempj=tempj+1;
                                    if newOD(i,tempj)-ODx_in_xy(i,tempj)<tempf-totempf
                                        totempf=totempf+newOD(i,tempj)-ODx_in_xy(i,tempj);
                                        xdown(1,j)=xdown(1,j)+newOD(i,tempj)-ODx_in_xy(i,tempj);
                                        ydown(1,j)=ydown(1,j)-(newOD(i,tempj)-ODx_in_xy(i,tempj));
                                        ODx_in_xy(i,tempj)=newOD(i,tempj);
                                    else
                                        ODx_in_xy(i,tempj)=ODx_in_xy(i,tempj)+(tempf-totempf);
                                        xdown(1,j)=xdown(1,j)+(tempf-totempf);
                                        ydown(1,j)=ydown(1,j)-(tempf-totempf);
                                        break;
                                    end
                                end
                                xup(1,i)=xup(1,i)+tempf;
                                yup(1,i)=yup(1,i)-tempf;
                            end
                        else
                            %%位置不够，有偏差
                            if persony+yup(1,i)<=max_py
                                %%大交路满，小交路有空
                                tempf=personx+xup(1,i)-max_px;
                                totempf=0;
                                personx=max_px;
                                persony=persony+yup(1,i);
                                deviaf=deviaf+tempf;   %累加记录真正溢出的偏差变量
                                if E==30
                                    break;
                                end
                                for j=E+1:30   %一定是坐到重合区间外的人坐不了车，否则他可以乘小交路
                                    px=floor(newOD(i,j)*tempf/xup(1,i));
                                    totempf=totempf+px;
                                    xdown(1,j)=xdown(1,j)-px;
                                    newOD(i,j)=newOD(i,j)-px;
                                end
                                tempj=E;
                                while totempf<tempf && tempj<30
                                    tempj=tempj+1;
                                    if newOD(i,tempj)<tempf-totempf
                                        totempf=totempf+newOD(i,tempj);
                                        xdown(1,j)=xdown(1,j)-newOD(i,tempj);
                                        newOD(i,tempj)=0;
                                    else
                                        newOD(i,tempj)=newOD(i,tempj)-(tempf-totempf);
                                        xdown(1,j)=xdown(1,j)-(tempf-totempf);
                                        break;
                                    end
                                end
                                xup(1,i)=xup(1,i)-tempf;
                            else
                                %%大小交路都满，此时偏差变量需要分别储存
                                tempfx=personx+xup(1,i)-max_px;
                                tempfy=persony+yup(1,i)-max_py;
                                totempfx=0;
                                totempfy=0;
                                personx=max_px;
                                persony=max_py;
                                deviaf=deviaf+tempfx+tempfy;
                                for j=i+1:30
                                    if j<=E
                                        py=floor(newOD(i,j)*tempfy/yup(1,i));
                                        totempfy=totempfy+py;
                                        ydown(1,j)=ydown(1,j)-py;
                                        newOD(i,j)=newOD(i,j)-py;
                                    else
                                        px=floor(newOD(i,j)*tempfx/xup(1,i));
                                        totempfx=totempfx+px;
                                        xdown(1,j)=xdown(1,j)-px;
                                        newOD(i,j)=newOD(i,j)-px;
                                    end
                                end
                                tempj=i;
                                while totempfy<tempfy && tempj<E
                                    tempj=tempj+1;
                                    if newOD(i,tempj)<tempfy-totempfy
                                        totempfy=totempfy+newOD(i,tempj);
                                        ydown(1,j)=ydown(1,j)-newOD(i,tempj);
                                        newOD(i,tempj)=0;
                                    else
                                        newOD(i,tempj)=newOD(i,tempj)-(tempfy-totempfy);
                                        ydown(1,j)=ydown(1,j)-(tempfy-totempfy);
                                        break;
                                    end
                                end
                                tempj=i;
                                while totempfx<tempfx && tempj<30
                                    tempj=tempj+1;
                                    if newOD(i,tempj)<tempfx-totempfx
                                        totempfx=totempfx+newOD(i,tempj);
                                        xdown(1,j)=xdown(1,j)-newOD(i,tempj);
                                        newOD(i,tempj)=0;
                                    else
                                        newOD(i,tempj)=newOD(i,tempj)-(tempfx-totempfx);
                                        xdown(1,j)=xdown(1,j)-(tempfx-totempfx);
                                        break;
                                    end
                                end
                                yup(1,i)=yup(1,i)-tempfy;
                                xup(1,i)=xup(1,i)-tempfx;
                            end
                        end
                    else   %不在重合区间
                        personx=personx-xdown(1,i);
                        if personx+xup(1,i)<=max_px
                            personx=personx+xup(1,i);
                        else
                            tempf=personx+xup(1,i)-max_px;
                            totempf=0;
                            personx=max_px;
                            deviaf=deviaf+tempf;
                            for j=i+1:30
                                px=floor(newOD(i,j)*tempf/xup(1,i));
                                totempf=totempf+px;
                                xdown(1,j)=xdown(1,j)-px;
                                newOD(i,j)=newOD(i,j)-px;
                            end
                            tempj=i;
                            while totempf<tempf && tempj<30
                                tempj=tempj+1;
                                if newOD(i,tempj)<tempf-totempf
                                    totempf=totempf+newOD(i,tempj);
                                    xdown(1,j)=xdown(1,j)-newOD(i,tempj);
                                    newOD(i,tempj)=0;
                                else
                                    newOD(i,tempj)=newOD(i,tempj)-(tempf-totempf);
                                    xdown(1,j)=xdown(1,j)-(tempf-totempf);
                                    break;
                                end
                            end
                            xup(1,i)=xup(1,i)-tempf;
                        end
                    end
                end

                deltaT=(xup/x+xdown/x+ydown/(x+y)+yup/(x+y))*0.04;   %停站时间,先根据总人数计算总停车时间
                
                interT=ceil(3600/(x+y));   %发车间隔
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
                deltaT(1,1)=0;
                deltaT(1,30)=0;
                
                proi=proi+1;   %方案满足约束，解数加一，并进行数据存储

                litkm=0;   %小交路公里数
                for i=B:E-1
                    litkm=litkm+dkm(i,1);
                end
                prokm(1,proi)=x*totalkm+y*litkm;   %总公里数

                totalT=0;   %总等待时间
                for i=1:29
                    for j=i+1:30
                        for di=i:j-1
                            totalT=totalT+newOD(i,j)*deltaT(1,di);
                        end
                        if i>=B && j<=E
                            totalT=totalT+newOD(i,j)*interT/2;
                        else
                            totalT=totalT+newOD(i,j)*interT*(x+y)/x/2;
                        end
                    end
                end
                proT(1,proi)=totalT/(sum(xup)+sum(yup));
                pros(1,proi)=x+y;
                prof(1,proi)=deviaf;

                prox(1,proi)=x;
                proE(1,proi)=E;
                proB(1,proi)=B;

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
                if prof(1,proi)<min_f
                    min_f=prof(1,proi);
                end
                if prof(1,proi)>max_f
                    max_f=prof(1,proi);
                end
            end
        end
    end
end
%%求变化范围
avers=(max_s-min_s);
averT=(max_T-min_T);
averkm=(max_km-min_km);
averf=(max_f-min_f);
%归一化平衡量纲
profit = (pros-min_s)/avers+(proT-min_T)*2/averT+(prokm-min_km)/averkm+(prof-min_f)/averf*4;   %使客流量的偏差变量占一半的比重
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
minx=prox(1,xi)   %19
miny=proy(1,xi)   %6
minE=proE(1,xi)   %17
minB=proB(1,xi)   %10
