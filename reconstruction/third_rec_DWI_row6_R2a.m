clc;clear all;close all;

%===========================
recsample(1,1)=struct('sample',[]);
load('./data/recsample_new/recsample62.mat');recsample(2,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample63.mat');recsample(3,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample64.mat');recsample(4,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample65.mat');recsample(5,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample66.mat');recsample(6,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample67.mat');recsample(7,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample68.mat');recsample(8,1)=struct('sample',sample);clear sample
load('./data/recsample_new/recsample69.mat');recsample(9,1)=struct('sample',sample);clear sample
load('../third_ADW_interpolation/data/china');
subchina=china_tele;clear china_tele china_notele
% china=china-3;

for ii=1:size(subchina,1)
    for jj=1:size(subchina,2)
        tt=squeeze(subchina(ii,jj,:));
        d1=mean(tt(1921-960+1:1980-960+1,1));
        d2=std(tt(1921-960+1:1980-960+1,1));
        
        for kk=1:length(tt)
            if ~isnan(subchina(ii,jj,kk))
                if tt(kk,1)<=(d1-1.17*d2)
                    subchina(ii,jj,kk)=1;
                elseif tt(kk,1)<=(d1-0.33*d2) & tt(kk,1)>(d1-1.17*d2)
                    subchina(ii,jj,kk)=2;
                elseif tt(kk,1)<=(d1+0.33*d2) & tt(kk,1)>(d1-0.33*d2)
                    subchina(ii,jj,kk)=3;
                elseif tt(kk,1)<=(d1+1.17*d2) & tt(kk,1)>(d1+0.33*d2)
                    subchina(ii,jj,kk)=4;
                elseif tt(kk,1)>(d1+1.17*d2);
                    subchina(ii,jj,kk)=5;
                end
                
            end
        end
        
        clear tt dd2 d1 d2
    end
end;clear ii jj kk
china=subchina;



load('./data/DWI_tree.mat');
load('../fourth_supplement_treering_in_misslocation/data/tree_row6.mat')
tree(6,4).tree=cat(2,tree(6,4).tree,row64);
tree(6,5).tree=cat(2,tree(6,5).tree,row65);

iii=6;
for ii=1:size(DWI,2)
    dd1=DWI(iii,ii).DWI;
    if ~isempty(dd1)
        for jj=1:size(dd1,2)
            dd2=dd1(:,jj);
            dd2_mean=mean(dd2(1921-960+1:1980-960+1));
            dd2_std=std(dd2(1921-960+1:1980-960+1));
            
            dd2=(dd2-dd2_mean)/dd2_std;
            dd1(:,jj)=dd2;
            clear dd2 dd2_mean dd2_std
        end
    end
    DWI(iii,ii).DWI=dd1;
    clear dd1
end;clear ii jj
        
for ii=1:size(tree,2)
    dd1=tree(iii,ii).tree;
    if ~isempty(dd1)
        for jj=1:size(dd1,2)
            dd2=dd1(:,jj);
            dd2_mean=mean(dd2(1921-960+1:1980-960+1));
            dd2_std=std(dd2(1921-960+1:1980-960+1));
            
            dd2=(dd2-dd2_mean)/dd2_std;
            dd1(:,jj)=dd2;
            clear dd2 dd2_mean dd2_std
        end
    end
    tree(iii,ii).tree=dd1;
    clear dd1
end;clear ii jj

val=[1951-960+1:1980-960+1]';
cal=[1921-960+1:1950-960+1]';
recon=[1921-960+1:1980-960+1]';

sss=6;
%============================reconstruction
for jj=1:size(china,2)
    if ~isempty(DWI(sss,jj).DWI) | ~isempty(tree(sss,jj).tree) 
        stationpre=cat(2,DWI(sss,jj).DWI,tree(sss,jj).tree);
        val_stationpre=stationpre(val,:);
        cal_stationpre=stationpre(cal,:);
        recon_stationpre=stationpre(recon,:);
        
        %=======================
        tic
        sample=recsample(jj,1).sample;
        for kk1=1:size(sample,1)
            for kk2=2:size(sample,2)
                kk33=size(sample(kk1,kk2).sample,2);
                for kk3=1:size(sample(kk1,kk2).sample,1)
                    %===================adjusted R2
                    xtree=recon_stationpre(:,sample(kk1,kk2).sample(kk3,:));
                    yDWI=squeeze(china(sss,jj,recon));
                    po=regress(yDWI,[ones(length(yDWI),1),xtree]);

                    yDWIX2=sum(repmat(po(2:end,1)',[size(xtree,1),1]).*xtree,2)+po(1,1);
                    [k,p]=corr(yDWI,yDWIX2);
                    r2(kk1,kk2).r2(kk3,1)=abs(((60-1)*k^2-size(xtree,2))/(60-size(xtree,2)-1));
                    
                    clear xtree yDWI po xtreeX yDWIX1 yDWIX2 k p
                    
                    %====================reconstructed DWI
                    xtree=stationpre(recon,sample(kk1,kk2).sample(kk3,:));
                    yDWI=squeeze(china(sss,jj,recon));
                    po=regress(yDWI,[ones(length(yDWI),1),xtree]);
                    
                    xtreeX=stationpre(sample(kk1,1).sample,sample(kk1,kk2).sample(kk3,:));
                    yDWIX2=sum(repmat(po(2:end,1)',[size(xtreeX,1),1]).*xtreeX,2)+po(1,1);
                    recdwi(kk1,kk2).sample(:,kk3)=yDWIX2;
                    
                    clear xtree yDWI po xtreeX yDWIX2
                    
                    %====================matching variance
                    xtree=stationpre(recon,sample(kk1,kk2).sample(kk3,:));
                    yDWI=squeeze(china(sss,jj,recon));
                    po=regress(yDWI,[ones(length(yDWI),1),xtree]);
                    
                    xtreeX=stationpre(recon,sample(kk1,kk2).sample(kk3,:));
                    yDWIX2=sum(repmat(po(2:end,1)',[size(xtreeX,1),1]).*xtreeX,2)+po(1,1);
                    match(kk1,kk2).sample(:,kk3)=yDWIX2;
                    
                    clear xtree yDWI po xtreeX yDWIX2
                    
                    %====================RE and CE: val:1951-1980
                    xtree=val_stationpre(:,sample(kk1,kk2).sample(kk3,:));
                    yDWI=squeeze(china(sss,jj,val));
                    po=regress(yDWI,[ones(length(yDWI),1),xtree]);
                    
                    xtreeX=stationpre(recon,sample(kk1,kk2).sample(kk3,:));
                    yDWIX1=squeeze(china(sss,jj,recon));
                    yDWIX2=sum(repmat(po(2:end,1)',[size(xtreeX,1),1]).*xtreeX,2)+po(1,1);
                    %----------RE
                    rece(kk1,kk2).rece(kk3,1)=1-sum((yDWIX1(31:end,1)-yDWIX2(31:end,1)).^2)/sum((yDWIX1(31:end,1)-mean(yDWIX1(31:end,1))).^2);
                    %----------CE
                    rece(kk1,kk2).rece(kk3,2)=1-sum((yDWIX1(1:30,1)-yDWIX2(1:30,1)).^2)/sum((yDWIX1(1:30,1)-mean(yDWIX1(1:30,1))).^2);

                    clear xtree yDWI po xtreeX yDWIX1 yDWIX2
                    
                    %====================RE and CE: cal:1951-1980
                    xtree=cal_stationpre(:,sample(kk1,kk2).sample(kk3,:));
                    yDWI=squeeze(china(sss,jj,cal));
                    po=regress(yDWI,[ones(length(yDWI),1),xtree]);
                    
                    xtreeX=stationpre(recon,sample(kk1,kk2).sample(kk3,:));
                    yDWIX1=squeeze(china(sss,jj,recon));
                    yDWIX2=sum(repmat(po(2:end,1)',[size(xtreeX,1),1]).*xtreeX,2)+po(1,1);
                    %----------RE
                    rece(kk1,kk2).rece(kk3,3)=1-sum((yDWIX1(1:30,1)-yDWIX2(1:30,1)).^2)/sum((yDWIX1(1:30,1)-mean(yDWIX1(1:30,1))).^2);
                    %----------CE
                    rece(kk1,kk2).rece(kk3,4)=1-sum((yDWIX1(31:end,1)-yDWIX2(31:end,1)).^2)/sum((yDWIX1(31:end,1)-mean(yDWIX1(31:end,1))).^2);

                    clear xtree yDWI po xtreeX yDWIX1 yDWIX2

                end
            end
        end;clear kk1 kk2 kk3
        recR2(1,jj)=struct('r2',r2);
        recDWI(1,jj)=struct('recDWI',recdwi);
        recRECE(1,jj)=struct('recRECE',rece);
        recMATCH(1,jj)=struct('recMATCH',match);
        clear r2 recdwi rece match
        
        disp(['The grid in row of ',num2str(sss),' and column of ',num2str(jj),' used ',num2str(toc),' seconds.'])
    end
    
    clear stationpre
end;clear jj


%===========================Best Subset Regression
for ii=1:size(china,2)
    if ~isempty(recR2(ii).r2)
        r2=recR2(ii).r2;
        for jj1=1:size(r2,1)
            for jj2=1:size(r2,2)
                if ~isempty(r2(jj1,jj2).r2)
                    dd1(jj2,1)=find(r2(jj1,jj2).r2==max(r2(jj1,jj2).r2));
                    dd1(jj2,2)=max(r2(jj1,jj2).r2);
                end
            end
            dd2=find(dd1(:,2)==max(dd1(:,2)));
            r2max(jj1,1)=dd2;
            r2max(jj1,2:3)=dd1(dd2,:);
                
            clear dd2 dd1
        end;clear jj1 jj2
            
        R2max(ii,1)=struct('rmax',r2max);
        clear r2max r2
    end
end;clear ii 

year=[1:size(china,3)]';
subchina=china(sss,:,:);
for ii=1:size(china,2)
    if ~isempty(R2max(ii).rmax)
        dd1=R2max(ii).rmax;
        r=cat(2,[1:size(dd1,1)]',dd1);
        subR2(1,ii)=mean(r(find(r(:,4)>0.05),4));
        rmaxnum=find(r(:,end)==max(r(:,end)));
        if length(rmaxnum)==1
            rmaxnum=rmaxnum;
        else
            rmaxnum=rmaxnum(end);
        end
        
        sample=recsample(ii).sample;
        DWI=recDWI(ii).recDWI;
        match=recMATCH(ii).recMATCH;
        
        rmax_mean=mean(match(r(rmaxnum,1),r(rmaxnum,2)).sample(:,r(rmaxnum,3)));
        rmax_std=std(match(r(rmaxnum,1),r(rmaxnum,2)).sample(:,r(rmaxnum,3)));
            
        %=========================match variance
        for jj=1:size(r,1)
            dd1=DWI(r(jj,1),r(jj,2)).sample(:,r(jj,3));
            ddd1=length(find(dd1-cat(1,dd1(2:end,1),dd1(end,1))==0));
            if jj~=rmaxnum & ddd1~=length(dd1)
                dd2(:,2)=(dd1-mean(dd1))/std(dd1)*rmax_std+rmax_mean;
            elseif jj==rmaxnum | ddd1==length(dd1)
                dd2(:,2)=dd1;
            end
            dd2(:,1)=sample(r(jj,1),1).sample;
            
            [z,ia,ib]=intersect(year,dd2(:,1));
            subchina(1,ii,ia)=dd2(:,2);
                
            clear dd1 dd2 dd3 z ia ib 
        end;clear jj
            
        clear dd1 r rmaxnum sample DWI match rmax_mean rmax_std rmaxnum
    end
end;clear ii


% subchina(subchina<-2)=-2;
% subchina(subchina>2)=2;
% subchina=round(subchina);
% subchina=subchina+3;

%============================transport into 5 grades
% for ii=1:size(subchina,1)
%     for jj=1:size(subchina,2)
%         tt=squeeze(subchina(ii,jj,:));
% %         if size(weight(ii,jj).CDD,1)~=1
%             tt1=zeros(size(tt,1),size(tt,2));
%             tt1(tt1==0)=nan;
% 
%             dd1=find(~isnan(tt));
%             dd2=tt(dd1,1);
%             dd2=zscore(dd2);
%             [d1,d2]=sort(dd2);
%             c1=floor(length(d2)*0.15);
%             c2=floor(length(d2)*0.35);
%             c3=floor(length(d2)*0.65);
%             c4=floor(length(d2)*0.85);
%     
%             tt1(dd1(d2(1:c1)),1)=1;
%             tt1(dd1(d2(c1+1:c2)),1)=2;
%             tt1(dd1(d2(c2+1:c3)),1)=3;
%             tt1(dd1(d2(c3+1:c4)),1)=4;
%             tt1(dd1(d2(c4+1:end)),1)=5;
%             
%             subchina(ii,jj,:)=tt1;
%             
%             clear tt1 dd1 dd2 c1 c2 c3 c4  
% %         end
%         
%         clear tt
%     end
% end;clear ii jj

for ii=1:size(subchina,1)
    for jj=1:size(subchina,2)
        tt=squeeze(subchina(ii,jj,:));
        d1=mean(tt(1921-960+1:1980-960+1,1));
        d2=std(tt(1921-960+1:1980-960+1,1));
        
        for kk=1:length(tt)
            if ~isnan(subchina(ii,jj,kk))
                if tt(kk,1)<=(d1-1.17*d2)
                    subchina(ii,jj,kk)=1;
                elseif tt(kk,1)<=(d1-0.33*d2) & tt(kk,1)>(d1-1.17*d2)
                    subchina(ii,jj,kk)=2;
                elseif tt(kk,1)<=(d1+0.33*d2) & tt(kk,1)>(d1-0.33*d2)
                    subchina(ii,jj,kk)=3;
                elseif tt(kk,1)<=(d1+1.17*d2) & tt(kk,1)>(d1+0.33*d2)
                    subchina(ii,jj,kk)=4;
                elseif tt(kk,1)>(d1+1.17*d2);
                    subchina(ii,jj,kk)=5;
                end
                
            end
        end
        
        clear tt dd2 d1 d2
    end
end;clear ii jj kk



for ii=1:size(subchina,2)
    dd1=squeeze(subchina(1,ii,:));
    submiss(1,ii)=length(find(isnan(dd1)))/size(subchina,3);
    
    clear dd1
end;clear ii
        
save('./data/sub_spe/zheng_subchina_row6_R2a.mat','subchina','subR2','recRECE','submiss','R2max');





















        
        
        
        
        
        


