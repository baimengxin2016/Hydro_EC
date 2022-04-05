clc;clear all;close all;

%============================
load('G:\IGSNRR\中国东部数据集重建\code\旱涝等级\重建\second_sensitivity\corr_EastChina\夏季\data\NHproxy_corr.mat')
sam=cat(1,sample(1).sample,sample(2).sample);
sam=cat(1,sam,sample(3).sample);
sam=sam(find(sam(:,4)<0),:);
uslonlat=sam(:,4:5);

%============================combine Northern Hemisphere TreeRing Index
load('G:\IGSNRR\中国东部数据集重建\proxydata\Asia\code\data\treering_Asia.mat');
asiatree=repmat(1,[length(960:2000) size(asiatree,1)]);
load('G:\IGSNRR\中国东部数据集重建\proxydata\Europe\code\data\treering_Euro.mat');
eurotree=repmat(1,[length(960:2000) size(eurotree,1)]);
load('G:\IGSNRR\中国东部数据集重建\code\全球代用资料特征\treering\corr_GPCC\data\treering_NA.mat');
natree=(natree(960:2000,:));
for ii=1:size(nayear,2)
    dd1=nayear(:,ii);
    dd2=find(~isnan(dd1));
    if length(dd2)~=0
        na(ii,1)=dd1(dd2(1),1);
    else
        na(ii,1)=nan;
    end
    
    clear dd1 dd2
end;clear ii
nayear=na;clear na

tree=cat(2,natree,asiatree);
tree=cat(2,tree,eurotree);
ustree=tree(:,sam(:,1));
rr=r(:,:,sam(:,1));
clear nalonlat natree nasampnum nayear asialonlat asiatree asiayear eurolonlat eurotree euroyear 


% %============================tranfer tree into DWI
% tt=ustree;
% for ii=1:size(ustree,2)
%     dd1=tt(1951-960+1:2000-960+1,ii);
%     d1=mean(dd1);
%     d2=std(dd1);
%     
%     for jj=1:size(tt,1)
%         if ~isnan(tt(jj,ii))
%             if tt(jj,ii)<=(d1-1.17*d2)
%                 tt1(jj,ii)=1;
%             elseif tt(jj,ii)<=(d1-0.33*d2) & tt(jj,ii)>(d1-1.17*d2)
%                 tt1(jj,ii)=2;
%             elseif tt(jj,ii)<=(d1+0.33*d2) & tt(jj,ii)>(d1-0.33*d2)
%                 tt1(jj,ii)=3;
%             elseif tt(jj,ii)<=(d1+1.17*d2) & tt(jj,ii)>(d1+0.33*d2)
%                 tt1(jj,ii)=4;
%             elseif tt(jj,ii)>(d1+1.17*d2);
%                 tt1(jj,ii)=5;
%             end
%             
%         else
%            tt1(jj,ii)=nan;
%            
%         end
%     end
%     clear d1 d2 dd1 dd2
% end;clear ii jj tt
% ustree=tt1;clear tt1

% for ii=1:size(ustree,2)
%     dd1=ustree(:,ii);
%     d1=mean(dd1(1951-960+1:2000-960+1,1));
%     d2=std(dd1(1951-960+1:2000-960+1,1));
%     
%     dd1=(dd1-d1)/d2;
%     ustree(:,ii)=dd1;
%     
%     clear dd1 d1 d2
% end;clear ii
    

for ii=1:size(rr,1)
    for jj=1:size(rr,2)
        dd1=squeeze(rr(ii,jj,:));
        dd2=find(~isnan(dd1));
        dd3=ustree(:,dd2);
        tr(ii,jj)=struct('tree',dd3);
        
        clear dd1 dd2 dd3
    end
end;clear ii jj
clear tree
tree=tr;clear tr

location=[1,1;2,1;4,1;5,1;6,1;7,1;8,1;
                 8,2;
                 8,7];
for ii=1:size(location,1)
    tree(location(ii,1),location(ii,2)).tree=[];
end;clear ii

for ii=1:size(tree,1)
    for jj=1:size(tree,2)
        dd1=tree(ii,jj).tree;
        if size(dd1,2)>=2
            dd1=dd1(:,find(~isnan(dd1(1,:))));
        end
        tree(ii,jj).tree=dd1;
        clear dd1
    end
end;clear ii jj


%============================DWI
fpath1=['G:\IGSNRR\中国东部数据集重建\code\旱涝等级\重建\third_ADW_interpolation\data\'];
load([fpath1,'DWI_historical_tree.mat'])

chinalonlat=lonlat;
load('G:\IGSNRR\中国东部数据集重建\code\旱涝等级\重建\third_ADW_interpolation\data\china.mat');clear china_notele

%============================find dry-wet index station 
for ii=1:length(lon)
    for jj=1:length(lat)
        x1=repmat(lon(ii,1),size(chinalonlat,1),1);
        y1=repmat(lat(jj,1),size(chinalonlat,1),1);
        
        %=========================the distance of grid point and station point
        x2=chinalonlat(:,1);y2=chinalonlat(:,2);
        distance=sqrt((x1-x2).^2+(y1-y2).^2);
        
        d1=find(distance<=400/111);
        d2=chinalonlat(d1,:);
        d2=cat(2,d2,d1);

        CDD(ii,jj)=struct('CDD',d2);
        
        clear x1 x2 distance d1 d2
    end
end;clear ii jj
for ii=1:size(CDD,1)
    for jj=1:size(CDD,2)
        dd1=DWI(:,CDD(ii,jj).CDD(:,3)+1);
        dw(ii,jj)=struct('DWI',dd1);
        clear dd1
    end
end;clear ii jj CDD
location=[1,1;2,1;4,1;5,1;6,1;7,1;8,1;
                 8,2;
                 8,7];
for ii=1:size(location,1)
    dw(location(ii,1),location(ii,2)).DWI=[];
end;clear ii
clear DWI
DWI=dw;


% save('./data/DWI_tree.mat','DWI','tree');




z=DWI(5,7).DWI;
zz=z(1:541,:);
year=[960:1500]';

for ii=1:size(zz,1)
    dd1=(zz(ii,:))';
    DD1=find(~isnan(dd1));
    DDD1=length(DD1);
    
    mm=0;
    for jj=ii+1:size(zz,1)
        dd2=(zz(jj,:))';
        DD2=find(~isnan(dd2));
        DDD2=length(DD2);
        if length(intersect(DD1,DD2))==max(DDD1,DDD2)
            mm=mm+1;
            zzz(ii,mm)=year(jj);
            continue
        end
        
        clear dd2 DD2
    end
    clear dd1 DD1
end;clear ii
zzz(:,1)=year;






