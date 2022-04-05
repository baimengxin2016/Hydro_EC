clc;clear all;close all;

%============================
load('./data/DWI_tree.mat');
% load('G:\IGSNRR\中国东部数据集重建\code\旱涝等级\重建\fourth_supplement_treering_in_misslocation\data\speleothems_row3.mat')
% tree(3,1).tree=cat(2,tree(3,1).tree,row31);

val=[1951-960+1:1980-960+1];
cal=[1921-960+1:1950-960+1];
%============================reconstruction
for ii=3%1:size(DWI,1)
    for jj=1:size(DWI,2)
        tic
        
        if ~isempty(DWI(ii,jj).DWI) | ~isempty(tree(ii,jj).tree) 
            stationpre=cat(2,DWI(ii,jj).DWI,tree(ii,jj).tree);
            val_stationpre=stationpre(val,:);
            cal_stationpre=stationpre(cal,:);
            
            if size(stationpre,2)>=2
                for kk=1:size(stationpre,1)
                    dd1=find(~isnan(stationpre(kk,:)));
                    if ismember([1:size(DWI(ii,jj).DWI,2)],dd1)
                        sam(kk,1)=struct('sam',[]);
                    else
                        sam(kk,1)=struct('sam',dd1);
                    end
                    clear dd1
                end;clear kk
                
                %=======================
                ss=sam;mm2=0;
                for kk1=1:length(ss)
                    if ~isempty(ss(kk1,1).sam)
                        mm1=0;
                        for kk2=1:length(ss)
                            if isequal(sam(kk2,1).sam,sam(kk1,1).sam)
                                mm1=mm1+1;
                                dd1(mm1,1)=kk2;
                                ss(kk2,1).sam=[];
                            end
                        end
                        mm2=mm2+1;
                        sample(mm2,1)=struct('sample',dd1);
                        clear dd1 mm1 
                    end
                end;clear kk1 kk2 ss mm2
                
                %=======================
                
                for kk1=1:length(sample)
                    site1=sample(kk1,1).sample(1,1);
                    site2=sam(site1,1).sam;
                    
                    for kk2=1:length(site2)
                        dd1=nchoosek(site2,kk2);
                        sample(kk1,kk2+1)=struct('sample',dd1);
                        
                        clear dd1
                    end;clear kk2
                    clear site1 site2
                end;clear kk1
                
                save(['./data/recsample_new/recsample',num2str(ii),num2str(jj),'.mat'],'sample','-v7.3');
                
                clear sample
            end
        end
        
        disp(['The grid in row of ',num2str(ii),' and column of ',num2str(jj),' used ',num2str(toc),' seconds.'])

    end
end;clear ii jj

























