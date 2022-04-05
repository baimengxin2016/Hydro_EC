clc;clear all;close all;

%============================historical document
fpath1=['G:\IGSNRR\中国东部数据集重建\code\旱涝等级\划分网格\data\'];
china=xlsread([fpath1,'DWI_LiuY.xlsx']);
lonlat=china(2:end,1:2);
DWI=china(:,3:end)';

%============================tree ring chronology width
fpath2=['G:\IGSNRR\中国东部数据集重建\code\旱涝等级\重建\first_treering_into_DWI\data\'];
load([fpath2,'treeDWI.mat']);
%============================combine historical document with tree ring
lonlat=cat(1,lonlat,treelonlat(:,1:2));
DWI=cat(2,DWI,treeDWI(:,2:end));
DWI(DWI==0)=nan;

% save('./data/DWI_historical_tree.mat','lonlat','DWI');

lon=[105:2.5:122.5]';
lat=[20:2.5:40]';

for ii=1:length(lon)
    for jj=1:length(lat)
        x1=repmat(lon(ii,1),size(lonlat,1),1);
        y1=repmat(lat(jj,1),size(lonlat,1),1);
        
        %=========================the distance of grid point and station point
        x2=lonlat(:,1);y2=lonlat(:,2);
        distance=sqrt((x1-x2).^2+(y1-y2).^2);
        
        d1=find(distance<=400/111);
        d2=lonlat(d1,:);
        d2=cat(2,d2,d1);

        CDD(ii,jj)=struct('CDD',d2);
        
        clear x1 x2 distance d1 d2
    end
end;clear ii jj

% %============================if grid only one station, the station must be in this grid
% for ii=1:size(CDD,1)
%     for jj=1:size(CDD,2)
%         dd1=CDD(ii,jj).CDD;
%         if size(dd1,1)==1
%             if abs(lon(ii,1)-dd1(1,1))>1.25 | abs(lat(jj,1)-dd1(1,2))>1.25
%                 dd1=[];
%             end
%         end
%         
%         CDD(ii,jj)=struct('CDD',dd1);
%         clear dd1
%     end
% end;clear ii jj

        
%============================interpolation
%=====first step: distance, named wk
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=CDD(ii,jj).CDD;
        x1=lon(ii,1);
        y1=lat(jj,1);
        
%         if ~isempty(dd1)
        if size(dd1,1)>=2
            for kk=1:size(dd1,1)
                x2=dd1(kk,1);
                y2=dd1(kk,2);
                dis=sqrt((x1-x2)^2+(y1-y2)^2);
                dd1(kk,4)=exp(-4*dis/(400/111));
                
                clear x2 y2 dis
            end
            
        end
        
        CDD(ii,jj)=struct('CDD',dd1);
        clear dd1 x1 y1
        
    end
end;clear ii jj kk

% stationnum=CDD;
% save('./data/stationnum.mat','stationnum');

%=====second step: directional (angular) isolation, named ak
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=CDD(ii,jj).CDD;
        x1=lon(ii,1);
        y1=lat(jj,1);
        
%         if ~isempty(dd1)
        if size(dd1)>=2
            for kk1=1:size(dd1,1)
                x2=dd1(kk1,1);y2=dd1(kk1,2);
                
                ss1=[1:size(dd1,1)]';
                ss2=ss1(ss1~=kk1);
              
                for kk2=1:length(ss2)
                    x3=dd1(ss2(kk2),1);y3=dd1(ss2(kk2),2);
                    
                    dis1=sqrt((x1-x2)^2+(y1-y2)^2);
                    dis2=sqrt((x1-x3)^2+(y1-y3)^2);
                    dis3=sqrt((x2-x3)^2+(y2-y3)^2);
                    cosA(kk2,1)=(dis1^2+dis2^2-dis3^2)/(2*dis1*dis2);
                    
                    clear  x3 y3 dis1 dis2 dis3
                end;
                
                ak=sum(dd1(ss2,4).*(1-cosA))/sum(dd1(ss2,4));
                dd1(kk1,5)=ak;
                
                clear x2 y2 ss1 ss2 cosA
            end;clear kk1 kk2
       
        end
        
        CDD(ii,jj)=struct('CDD',dd1);
        clear dd1
    end
end;clear ii jj 

%=====third step: the angular and distance weights
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=CDD(ii,jj).CDD;
%         if ~isempty(dd1)
        if size(dd1,1)>=2
            dd1(:,6)=dd1(:,4).*(1+dd1(:,5));
            CDD(ii,jj)=struct('CDD',dd1);
        
            clear dd1
        end
    end
end;clear ii jj

%============================
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=CDD(ii,jj).CDD;
        if size(dd1,1)>=2
            dd1(:,6)=dd1(:,6)/sum(dd1(:,6));
            weight(ii,jj)=struct('CDD',dd1);
        
            clear dd1
        end
    end
end;clear ii jj

%============================if grid only one station, the station is the grid
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=CDD(ii,jj).CDD;
        if size(dd1,1)==1
            weight(ii,jj)=struct('CDD',cat(2,CDD(ii,jj).CDD,[1,1,1]));
        end
        
        clear dd1
    end
end;clear ii jj


% save('./data/weight.mat','weight');


for ii=1:size(weight,1)
    for jj=1:size(weight,2)
        if ~isempty(weight(ii,jj).CDD)
            z(ii,jj)=sum(weight(ii,jj).CDD(:,6));
        else
            z(ii,jj)=nan;
        end
    end
end;clear ii jj




for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        station(ii,jj)=size(CDD(ii,jj).CDD,1);
    end
end;clear ii jj











