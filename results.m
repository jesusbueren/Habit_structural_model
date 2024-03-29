clear all
% close all
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')
fileID = fopen('asset_distribution_model.txt','r');
A=fscanf(fileID,'%f');
fclose(fileID);
A = standardizeMissing(A,-9)
T=36
type_y=2
type_e=3
model=reshape(A,7,T,type_y,type_e)

fileID = fopen('delta_assets_model_vs_data.txt','r');
A=fscanf(fileID,'%f');
fclose(fileID);
A = standardizeMissing(A,-9)
T=36
delta_a=reshape(A,6,T,type_e)

cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')

fileID=fopen('wealth_moments_data.txt');
mean_wealth=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);

data=reshape(mean_wealth{1},9,37,2,type_y,3);

colors = {  [0,   0,   0.2]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
colors3 =   [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .4 .4 .4] ;
lw=[1.5 1.5 2.5]
FS=10
sz=30
h_l=1
ini_age=26
fin_age=76 
last=size(ini_age:2:fin_age,2)

figure(1)
set(1,'position',[150    150    750    400]) %get(0, 'Screensize')
for y_l=1:type_y
for e_l=1:type_e
    h((e_l-1)*3+y_l)=subplot(3,2,(e_l-1)*2+y_l)
    C=plot(ini_age:2:fin_age,model(4,1:last,y_l,e_l),'Color',colors{1},'linewidth',1.5,'linestyle',pattern{1})
    hold on
    A=plot(ini_age:2:fin_age,model(5,1:last,y_l,e_l),'Color',colors{1},'linewidth',1.5,'linestyle',pattern{2},'HandleVisibility','off')
    plot(ini_age:2:fin_age,model(6,1:last,y_l,e_l,h_l),'Color',colors{1},'linewidth',1.5,'linestyle',pattern{1},'HandleVisibility','off')
    D=scatter(ini_age:4:fin_age,data(5,1:2:last,h_l,y_l,e_l)./1000,sz,colors3(3,:),'filled','s')
    B=scatter(ini_age:4:fin_age,data(6,1:2:last,h_l,y_l,e_l)./1000,sz,colors3(3,:),'filled','d','HandleVisibility','off')
    scatter(ini_age:4:fin_age,data(7,1:2:last,h_l,y_l,e_l)./1000,sz,colors3(3,:),'filled','s','HandleVisibility','off')
    if e_l==1 && y_l==1
        title('Protective')
%         ylim([0 150])
    elseif e_l==1 && y_l==2
        title('Detrimental')
%         ylim([0 600])
    elseif e_l==1 && y_l==3
        title('Harmful')
%         ylim([0 1000])
    end
    if y_l==1
        if e_l==1
            ylabel('Dropout')
        elseif e_l==2
            ylabel('High-school')
        elseif e_l==3
            ylabel('College')
        end 
    end
    if e_l==1
        ylim([0 700])
        yticks(0:200:700)
    elseif e_l==2
        ylim([0 1000])
        yticks(0:300:1000)
    else
        ylim([0 1800])
        yticks(0:500:1800)
    end
%     if e_l==1
%         ylim([0 200])
%         yticks(0:50:200)
%     elseif e_l==2
%         ylim([0 200])
%         yticks(0:50:200)
%     else
%         ylim([0 200])
%         yticks(0:50:200)
%     end
    xlim([ini_age-1 fin_age+1])
    set(gca,'FontName','Times New Roman','FontSize',FS);
    set(gcf,'color','w')
end
end

% I=legend([A B C D],'p50: Model','p50: Data','p25-p75: Model','p25-p75: Data','Location','northwest','orientation','horizontal')
% legend('boxoff')
% I.FontSize=FS
% newPosition = [0.45 -0.02 0.1 0.1];
% newUnits = 'normalized';
% set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\model_fit','-depsc')

% for y_l=1:type_y
% for e_l=1:type_e
%     if ((e_l-1)*3+y_l>1)
%         close(2)
%     end
%     hfig = figure(2);
%     set(2,'position',[150    150    380    200])
%     hax_new = copyobj(h((e_l-1)*3+y_l), hfig);
%     set(gcf,'color','w')
%     if e_l==1
%         edu_title="Dropout"
%     elseif e_l==2
%         edu_title="Highschool"
%     elseif e_l==3
%         edu_title="College"
%     end 
%     if y_l==1
%         type_title="Protective"
%     elseif y_l==2
%         type_title="Detrimental"
%     elseif y_l==3
%         type_title="Harmful"
%     end 
%     ylabel('$000')
%     legend('model','data','Location','northwest','Box','off','FontSize',FS+2)
%     set(gca,'FontName','Times New Roman','FontSize',FS);
%     str_title=strcat(edu_title," ", type_title)
%     str_title2=strcat(edu_title, type_title)
%     title(str_title)
%     set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
%     print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\model_fit',str_title2),'-depsc')
%     set(gca,'FontName','Times New Roman','FontSize',FS);
% end
% end



% figure(3)
% set(3,'position',[150    150    500    400])
% for e_l=1:3
%     for h_l=1:2
%         subplot(3,2,h_l+(e_l-1)*2)
%         plot(52:2:90,delta_a(2+h_l,14:33,e_l),'Color',colors{1},'linewidth',1.5,'linestyle',pattern{2})
%         hold on
%         scatter(52:2:90,delta_a(4+h_l,14:33,e_l),sz,colors3(3,:),'filled','d','HandleVisibility','off')
%         ylim([-50 50])
%         if h_l==1
%             if e_l==1
%                 ylabel('Dropout')
%             elseif e_l==2
%                 ylabel('High-school')
%             elseif e_l==3
%                 ylabel('College')
%             end 
%         end
%         if e_l==1
%             if h_l==1
%                 title('Good health')
%             elseif h_l==2
%                 title('Bad health')
%             end 
%         end
%         set(gca,'FontName','Times New Roman','FontSize',FS);
%         set(gcf,'color','w')
%     end
% end

% print('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\model_fit2','-depsc')





%% Costs
clear all
close all
clc
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')
data=importdata('costs_e_y.txt')

colors = {  [0.4660    0.6740    0.1880]+0.1   [0.9290    0.6940    0.1250]+0.1    [0.8500    0.3250    0.0980]+0.1 };
colors2 ={  [0.4660    0.6740    0.1880]-0.2   [0.9290    0.6940    0.1250]-0.2    [0.8500    0.3250    0.0980]-0.2   };
colors3 ={  [0.4660    0.6740    0.1880]-0.3   [0.9290    0.6940    0.1250]-0.2    [0.8500    0.3250    0.0980]-0.3   };
colors_s =   [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250] ; %; 0.8500    0.3250    0.0980
FS=10
type_y=2
G_df=1
data_mat=data.data
data_mat=reshape(data_mat(:,2:5),type_y,3,G_df+1,4)
% counterfactual=counterfactual.data
% counterfactual=reshape(counterfactual(:,2:4),3,3,2,3)
costs=data_mat(:,:,:,1)


for j=1:2
    figure(j)
    set(j,'position',[150    150    500    300])
for df_l=1:G_df
for e_l=1:3


X=categorical({'Protective','Detrimental'})
X = reordercats(X,{'Protective','Detrimental'});

if j==1 
    subplot(G_df,3,e_l+(df_l-1)*3)
    b=bar(X,data_mat(:,e_l,df_l,j),'FaceColor','flat')
%     for c_l=1:3
%         b.CData(c_l,:) = colors{c_l}
%     end
elseif j==2 
    subplot(G_df,3,e_l+(df_l-1)*3)
    b=bar(X,data_mat(:,e_l,df_l,j),'FaceColor','flat')
%     for c_l=1:3
%         b.CData(c_l,:) = colors{c_l}
%     end
end
hold on
if e_l==1
    title('Dropouts')
elseif e_l==2 
    title('Highschool')
elseif e_l==3
    title('College')
end
if e_l==1
    if df_l==1
        ylabel('Impatient')
    elseif df_l==2
        ylabel('Patient')
    elseif df_l==3
        ylabel('All')
    end
end

set(gca,'FontName','Times New Roman','FontSize',FS);
if j==1
    ylim([0,30])
elseif j==2
   ylim([min(data_mat(:,:,df_l,j),[],'all')-10,max(data_mat(:,:,df_l,j),[],'all')+10])

elseif j==3
    ylim([0,1])
end
set(gcf,'color','w')
end
print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\stage1_',num2str(j)),'-depsc')
end
end


for df_l=1:G_df
 j=3

X=categorical({'Dropouts','Highschool','College'})
X = reordercats(X,{'Dropouts','Highschool','College'});

figure(3)
set(3,'position',[150    150    700    300])


subplot(1,G_df,df_l)
for e_l=1:3
    A(e_l,:,:)=squeeze(data_mat(:,e_l,df_l,j:j+1))./sum(squeeze(data_mat(:,e_l,df_l,j:j+1)))
    B(e_l,:)=A(e_l,2,:)
end
b=bar(X,B,'FaceColor','flat') 
%     for c_l=1:3
%        b(1).CData(c_l,:) = colors{c_l}
%        b(2).CData(c_l,:) = colors2{c_l}
%        b(3).CData(c_l,:) = colors3{c_l}
%     end

% if df_l==1
%     title('Impatient')
% elseif df_l==2 
%     title('Patient')
% elseif df_l==3
%     title('All')
% end
legend('Data','Model')
set(gca,'FontName','Times New Roman','FontSize',FS);
if j==1
%     ylim([-1.5,1.5])
elseif j==2
   ylim([min(data_mat(:,:,df_l,j),[],'all')-10,max(data_mat(:,:,df_l,j),[],'all')+10])

elseif j==3
    ylim([0,0.7])
end
set(gcf,'color','w')
end
% print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\stage1_',num2str(j)),'-depsc')




for df_l=1:2
    df_l
data_mat(1,1,df_l,2)-data_mat(2,1,df_l,2)
data_mat(1,3,df_l,2)-data_mat(2,3,df_l,2)
data_mat(1,3,df_l,2)-data_mat(2,3,df_l,2)-(data_mat(1,1,df_l,2)-data_mat(2,1,df_l,2))
end


%%
sz = 50;
for df_l=1:2
    x=data_mat(:,:,df_l,2)
    y=log(data_mat(:,:,df_l,3))
    figure(4)
    set(4,'position',[150    150    450    200])
    subplot(1,2,df_l)
    h1=scatter(x(:),y(:))
    hl=lsline
    hl.LineWidth=2
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    Slope = B(2)
    Intercept = B(1)
    e_hat=y(:)-(Intercept+Slope.*x(:))
    Rsq=1-sum(e_hat.^2)/sum((y(:)-mean(y(:))).^2)
    delete(h1)
    hold on
    s=scatter(x(:,1),y(:,1),sz,colors_s,'filled','d')
    s=scatter(x(:,2),y(:,2),sz+10,colors_s,'filled','s')
    s=scatter(x(:,3),y(:,3),sz,colors_s,'filled','o')
    qw{1} = scatter(nan,nan,sz,[0 0 0], 'filled','d');
    qw{2} = scatter(nan,nan,sz+10,[0 0 0], 'filled','s');
    qw{3} = scatter(nan,nan,sz,[0 0 0], 'filled','o');
    legend([qw{1},qw{2},qw{3}], {'DO','HS','Co'}, 'location', 'best','Orientation','horizontal','Box','off','Position',[0.4 0.5 0.2 0.96])
    xlabel('Lifetime Utility')
    ylabel('log share')

    ylim([-7 0])

%     xlim([5 16])
    theString = sprintf('log(share) = %.3f V_0 %.2f ', Slope, Intercept);
    xL=xlim
    yL=ylim
    text((xL(2)-xL(1))*0.1+xL(1), (yL(2)-yL(1))*0.17+yL(1), theString,'FontName','Times New Roman', 'FontSize', FS)
    theString = sprintf('R^2= %.2f',Rsq);
    text((xL(2)-xL(1))*0.1+xL(1), (yL(2)-yL(1))*0.07+yL(1), theString,'FontName','Times New Roman', 'FontSize', FS)
    set(gca,'FontName','Times New Roman','FontSize',FS);
    set(gcf,'color','w')
end
print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\stage1'),'-depsc')



