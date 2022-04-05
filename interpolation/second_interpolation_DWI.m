clc;clear all;close all;

%============================
load('./data/DWI_historical_tree.mat')
DWI(:,2:end)=DWI(:,2:end)-3;%---anomaly

load('./data/weight.mat');

lon=[105:2.5:122.5]';
lat=[20:2.5:40]';
year=[960:2000]';

china1=zeros(length(lon),length(lat),length(year));
china1(china1==0)=nan;
china2=zeros(length(lon),length(lat),length(year));
china2(china2==0)=nan;

for ii=1:size(weight,1)
    for jj=1:size(weight,2)
        dd1=weight(ii,jj).CDD;
        
        if size(dd1,1)>=2
            china1(ii,jj,:)=sum(DWI(:,dd1(:,3)+1).*repmat((dd1(:,6))',size(china1,3),1),2);
            china2(ii,jj,:)=nansum(DWI(:,dd1(:,3)+1).*repmat((dd1(:,6))',size(china2,3),1),2);
            
        elseif size(dd1,1)==1
            china1(ii,jj,:)=DWI(:,dd1(1,3));
            china2(ii,jj,:)=DWI(:,dd1(1,3));
            
        end
        
        clear dd1
    end
end;clear ii jj

% china1(china1<-2)=-2;
% china1(china1>2)=2;
% china1=round(china1);
% china1=china1+3;
% 
% china2(china2<-2)=-2;
% china2(china2>2)=2;
% china2=round(china2);
% china2=china2+3;

% for ii=1:size(china1,1)
%     for jj=1:size(china1,2)
%         tt=squeeze(china1(ii,jj,:));
%         if size(weight(ii,jj).CDD,1)~=1
%             
%             d1=mean(tt(1951-960+1:2000-960+1,1));
%             d2=std(tt(1951-960+1:2000-960+1,1));
%         
%             tt1=zeros(size(tt,1),1);
%             tt1(tt1==0)=nan;
%             for kk=1:size(tt,1)
%                 if tt(kk,1)<=(d1-1.17*d2)
%                     tt1(kk,1)=1;
%                 elseif tt(kk,1)<=(d1-0.33*d2) & tt(kk,1)>(d1-1.17*d2)
%                     tt1(kk,1)=2;
%                 elseif tt(kk,1)<=(d1+0.33*d2) & tt(kk,1)>(d1-0.33*d2)
%                     tt1(kk,1)=3;
%                 elseif tt(kk,1)<=(d1+1.17*d2) & tt(kk,1)>(d1+0.33*d2)
%                     tt1(kk,1)=4;
%                 elseif tt(kk,1)>(d1+1.17*d2);
%                     tt1(kk,1)=5;
%                 end
%             end;clear kk
%         
%             china1(ii,jj,:)=tt1;
%             clear tt tt1 d1 d2
%             
%         end
%     end
% end;clear ii jj


location=[1,1;2,1;4,1;5,1;6,1;7,1;8,1;
                 8,2;
                 8,7];
for ii=1:size(location,1)
    china1(location(ii,1),location(ii,2),:)=nan;
end;clear ii


for ii=1:size(location,1)
    china2(location(ii,1),location(ii,2),:)=nan;
end;clear ii


china_tele=china1;
china_notele=china2;

% save('./data/china.mat','china_tele','china_notele','lon','lat');


for ii=1:size(china_tele,1)
    for jj=1:size(china_tele,2)
        figure;
        dd1=squeeze(china_tele(ii,jj,:));
        plot([960:2000]',dd1);
        
    end
end;clear ii jj

        
        
clc;clear all;close all;

load('./data/china.mat');
china=china_tele;

for ii=1:8
    zz=squeeze(china(ii,:,:));
    figure;
    
    for jj=1:9
        subplot(3,3,jj);
        plot([960:2000]',squeeze(zz(jj,:)))
        
        set(gca,'xlim',[950 2010],'xtick',[1000:200:2000])
        
    end
        
end;clear ii jj

        
        
        
        
        
        
        
        
