clear all
close all
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
%     if e_l==1
%         ylim([0 700])
%         yticks(0:200:700)
%     elseif e_l==2
%         ylim([0 1000])
%         yticks(0:300:1000)
%     else
        ylim([0 1800])
        yticks(0:500:1800)
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

for y_l=1:type_y
for e_l=1:type_e
    if ((e_l-1)*3+y_l>1)
        close(2)
    end
    hfig = figure(2);
    set(2,'position',[150    150    380    200])
    hax_new = copyobj(h((e_l-1)*3+y_l), hfig);
    set(gcf,'color','w')
    if e_l==1
        edu_title="Dropout"
    elseif e_l==2
        edu_title="Highschool"
    elseif e_l==3
        edu_title="College"
    end 
    if y_l==1
        type_title="Protective"
    elseif y_l==2
        type_title="Detrimental"
    elseif y_l==3
        type_title="Harmful"
    end 
    ylabel('$000')
    legend('model','data','Location','northwest','Box','off','FontSize',FS+2)
    set(gca,'FontName','Times New Roman','FontSize',FS);
    str_title=strcat(edu_title," ", type_title)
    str_title2=strcat(edu_title, type_title)
    title(str_title)
    set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
    print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\model_fit',str_title2),'-depsc')
    set(gca,'FontName','Times New Roman','FontSize',FS);
end
end



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





%% First-stage: model fit
clear all
close all
% clc
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')
data=importdata('costs_e_y.txt')

colors = [ 0.0000, 0.4470, 0.7410;   % Dark Blue
           0.8500, 0.3250, 0.0980;   % Dark Red
           0.4660, 0.6740, 0.1880;   % Dark Green
           0.3, 0.1, 0.4];  % Dark Purple
FS=10;
type_y=2;
G_df=1;
G_cohorts=2;
data_mat=data.data;
data_mat=reshape(data_mat(:,2:5),type_y,3,G_cohorts,G_df,4);

e_l=1
c_l=1
df_l=1
j=2
data_mat(1,e_l,c_l,df_l,j)-data_mat(2,e_l,c_l,df_l,j)
e_l=1
c_l=2
df_l=1
j=2
data_mat(1,e_l,c_l,df_l,j)-data_mat(2,e_l,c_l,df_l,j)

% (data_mat(1,1,1,1,2)-data_mat(2,1,1,1,2))-((data_mat(1,1,2,1,2)-data_mat(2,1,2,1,2)))
% 
% (data_mat(1,3,1,1,2)-data_mat(2,3,1,1,2))-((data_mat(1,3,2,1,2)-data_mat(2,3,2,1,2)))



(data_mat(1,3,1,1,2)-data_mat(1,1,1,1,2))/(data_mat(1,2,1,1,2)-data_mat(1,1,1,1,2))

% counterfactual=counterfactual.data
% counterfactual=reshape(counterfactual(:,2:4),3,3,2,3)
costs=data_mat(:,:,:,1)

c_l=1

for j=1:3
    figure(j)
    set(j,'position',[150    150    500    300])
for df_l=1:G_df
for e_l=1:3

X=categorical({'Protective','Detrimental'})
X = reordercats(X,{'Protective','Detrimental'});

if j==1 
    subplot(G_df,3,e_l+(df_l-1)*3)
    b=bar(X,data_mat(:,e_l,c_l,df_l,j),'FaceColor','flat')
%     for c_l=1:3
%         b.CData(c_l,:) = colors{c_l}
%     end
elseif j==2 
    subplot(G_df,3,e_l+(df_l-1)*3)
    b=bar(X,[data_mat(:,e_l,c_l,df_l,j) data_mat(:,e_l,2,df_l,j)],'FaceColor','flat')
%     for c_l=1:3
%         b.CData(c_l,:) = colors{c_l}
%     end
elseif j==3 
    subplot(G_df,3,e_l+(df_l-1)*3)
    b=bar(X,[data_mat(:,e_l,c_l,df_l,2)-data_mat(:,e_l,c_l,df_l,1) data_mat(:,e_l,2,df_l,2)-data_mat(:,e_l,c_l,df_l,1)],'FaceColor','flat')
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
    ylim([-5,100])
elseif j==2
   ylim([min(data_mat(:,:,c_l,df_l,j),[],'all')-10,max(data_mat(:,:,c_l,df_l,j),[],'all')+10])

end
set(gcf,'color','w')
end
print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\2023_EUI\figures\stage1_',num2str(j)),'-depsc')
end
end


df_l=1
for c_l=1:G_cohorts
j=3
X=categorical({'Dropouts','Highschool','College'})
X = reordercats(X,{'Dropouts','Highschool','College'});
figure(4)
set(4,'position',[150    150    500    250])
subplot(1,G_cohorts,c_l)
for e_l=1:3
        A(e_l,:,:)=squeeze(data_mat(:,e_l,c_l,df_l,j:j+1))./sum(squeeze(data_mat(:,e_l,c_l,df_l,j:j+1)));
        B(e_l,:)=A(e_l,2,:);
end
b=bar(X,B,'FaceColor','flat')
if c_l==1
    legend('Data','Model','Location','northwest')
end
set(gca,'FontName','Times New Roman','FontSize',FS-1);
if c_l==1
    title('Born in 1930')
elseif c_l==2
    title('Born in 1970')
elseif c_l==3
    title('Born in 1950')
elseif c_l==4
    title('Born in 1970')
end
ylim([0,1])
set(gcf,'color','w')
end
print('C:\Users\jbueren\Dropbox\habits\Slides\2023_EUI\figures\stage1_pr_y_e','-depsc')


df_l=1;

for c_l=1:G_cohorts
    j=3;
    X={'Dropouts','Highschool','College','Protective','Detrimental'};

    figure(5);
    set(5,'position',[150 150 500 250]);
    subplot(1,G_cohorts,c_l);

    for e_l=1:3
        A(e_l,:,:)=squeeze(data_mat(:,e_l,c_l,df_l,j:j+1));
        B2(e_l,:)=sum(squeeze(A(e_l,:,:)),1);
    end
    for e_l=1:2
        B2(e_l+3,:)=sum(squeeze(data_mat(e_l,:,c_l,df_l,j:j+1)));
    end

    % Use numeric values as x-data for plotting bars
    xVals = [1, 2, 3, 4.5, 5.5]; % Introduce a gap after 'College'
    
    for i = 1:numel(X)
        b1=bar(xVals(i), B2(i, :), 'FaceColor', 'flat');
        hold on;
        set(b1(1), 'FaceColor', colors(1, :));
        set(b1(2), 'FaceColor', colors(2, :));
    end
    hold off;
    % Adjusting x-ticks and labels
    set(gca, 'XTick', xVals);
    set(gca, 'XTickLabel', X);
    xline(3.75,'--'); % Adjust the xline to account for the modified x-values
    if c_l==1
        legend([b1(1),b1(2)],'Data','Model','Location','northwest')
    end
    set(gca,'FontName','Times New Roman','FontSize',FS-1);
    ylim([0,1.0]);
    if c_l==1
        title('Born in 1930');
    elseif c_l==2
        title('Born in 1970');
    elseif c_l==3
        title('Born in 1950');
    elseif c_l==4
        title('Born in 1970');
    end
    set(gcf,'color','w');
end

print('C:\Users\jbueren\Dropbox\habits\Slides\2023_EUI\figures\stage1_pr_e','-depsc');







% 
% df_l=1
% for c_l=1:G_cohorts
% j=3
% X=categorical({'Dropouts','Highschool','College','Protective','Detrimental'})
% X = reordercats(X,{'Dropouts','Highschool','College','Protective','Detrimental'});
% figure(4)
% set(4,'position',[150    150    500    250])
% subplot(1,G_cohorts,c_l)
% for e_l=1:3
%     A(e_l,:,:)=squeeze(data_mat(:,e_l,c_l,df_l,j:j+1))
%     B2(e_l,:)=sum(squeeze(A(e_l,:,:)),1)
% end
% for e_l=1:2
%     B2(e_l+3,:)=sum(squeeze(data_mat(e_l,:,c_l,df_l,j:j+1)))
% end
% b=bar(X,B2,'FaceColor','flat') 
% 
% set(gca,'FontName','Times New Roman','FontSize',FS-1);
% xline(3.5,'--');
% if c_l==1
% legend(b,'Data','Model','Location','northwest')
% end
% ylim([0,1.0])
% if c_l==1
%     title('Born in 1930')
% elseif c_l==2
%     title('Born in 1950')
% elseif c_l==3
%     title('Born in 1970')
% end
% set(gcf,'color','w')
% end
% print('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\stage1_pr_e','-depsc')

%% Decomposition first stage
clear all
close all
clc
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')

colors = [ 0.0000, 0.4470, 0.7410;   % Dark Blue
           0.8500, 0.3250, 0.0980;   % Dark Red
           0.4660, 0.6740, 0.1880;   % Dark Green
           0.3, 0.1, 0.4];  % Dark Purple




FS=10
type_y=2
G_df=1
benchmark=importdata('costs_e_y.txt')
benchmark=benchmark.data
benchmark=reshape(benchmark(:,2:5),type_y,3,2,G_df,4)


no_corr=importdata('costs_e_y_no_corr.txt')
no_corr=no_corr.data
no_corr=reshape(no_corr(:,2:5),type_y,3,2,G_df,4)

c_H=importdata('costs_e_y_H.txt')
c_H=c_H.data
c_H=reshape(c_H(:,2:5),type_y,3,2,G_df,4)

c_I=importdata('costs_e_y_I.txt')
c_I=c_I.data
c_I=reshape(c_I(:,2:5),type_y,3,2,G_df,4)


df_l=1
for c_l=[1 2]
j=3
X=categorical({'Dropouts','College'})
X = reordercats(X,{'Dropouts','College'})
for e_l=[1 3]

    A(e_l,:,:)=squeeze([benchmark(:,e_l,c_l,df_l,j+1) no_corr(:,e_l,c_l,df_l,j+1) c_I(:,e_l,c_l,df_l,j+1) c_H(:,e_l,c_l,df_l,j+1)  ])./sum(squeeze([benchmark(:,e_l,c_l,df_l,j+1) no_corr(:,e_l,c_l,df_l,j+1) c_I(:,e_l,c_l,df_l,j+1) c_H(:,e_l,c_l,df_l,j+1) ]));
    B(e_l,:)=A(e_l,2,:);

end
figure(c_l)
set(c_l,'position',[150    150    500    250])
subplot(1,2,2)
% Plot the Dropouts group
b1 = bar(X(1), B(1,:));
hold on
% Plot the College group
b2 = bar(X(2), B(3,:));

% Assuming colors is your RGB matrix with a color for each bar.
for k = 1:length(b1)
    set(b1(k), 'FaceColor', colors(k, :));
    set(b2(k), 'FaceColor', colors(k, :));
end

legend('Benchark','No correlation (\epsilon_e,\epsilon_y)', 'Equal Income', 'Equal Health' )

set(gca,'FontName','Times New Roman','FontSize',FS-2);
% if c_l==1
%     title('Born in 1930')
% else
%     title('Born in 1970')
% end
ylim([0,1])
set(gcf,'color','w')
%end


clear A B2
df_l=1
%for c_l=[1 3]
j=3
X=categorical({'Dropouts','Highschool','College','Protective','Detrimental'})
X = reordercats(X,{'Dropouts','Highschool','College','Protective','Detrimental'});

for e_l=1:3
    A(e_l,:,:)=squeeze([benchmark(:,e_l,c_l,df_l,j+1) no_corr(:,e_l,c_l,df_l,j+1) c_I(:,e_l,c_l,df_l,j+1) c_H(:,e_l,c_l,df_l,j+1)])
    B2(e_l,:)=sum(squeeze(A(e_l,:,:)),1)
end
for e_l=1:2
    B2(e_l+3,:)=sum(squeeze([benchmark(e_l,:,c_l,df_l,j+1); no_corr(e_l,:,c_l,df_l,j+1); c_I(e_l,:,c_l,df_l,j+1); c_H(e_l,:,c_l,df_l,j+1)]),2)
end
figure(c_l)
set(c_l,'position',[150    150    500    250])
subplot(1,2,1)
b = bar(X,B2,'FaceColor','flat'); 
    for k = 1:size(B2, 2)  % loop over the groups of bars
        b(k).CData = repmat(colors(k,:), [size(B2, 1), 1]);  % Apply color to each bar in the group
    end

yticks(0:0.2:1)
set(gca,'FontName','Times New Roman','FontSize',FS-2);
xline(3.5,'--');

ylim([0,1])
% if c_l==1
%     title('Born in 1930')
% else
%     title('Born in 1970')
% end
set(gcf,'color','w')
print(strcat('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\counterfactual_uncond',string(c_l)),'-depsc')
end


%% distribution of epsilon_pro across educations
clear all
close all
clc
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')

colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};



FS=10
G_educ=3
benchmark=importdata('eps_HSD.txt');
[f,xi]=ksdensity(benchmark(:,3));
figure(1)
set(1,'position',[150    150    400    250])
plot(xi,f,'Linewidth',2)
xline(min(benchmark(benchmark(:,1)==1,3)),'-','Linewidth',2)
xline(min(benchmark(benchmark(:,2)==1,3)),'--','Linewidth',2)
legend('pdf ','1930','1970')
xlim([-1.5 1.5])
set(gca,'FontName','Times New Roman','FontSize',FS-2);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\2023_EUI\figures\eps_pro_hsd','-depsc');


%%

c_l=1
figure(1)

   subplot(1,2,min(e_l,2))
   for y_l=1:2
    
    if y_l==1
       xline(min(benchmark((benchmark(2,:,c_l)==e_l &  benchmark(3,:,c_l)==y_l),c_l)),'-','1930','Linewidth',2)
       xline(min(benchmark(1,(benchmark(2,:,2)==e_l &  benchmark(3,:,2)==y_l),2)),'-','1970','Linewidth',2)
    end
    if e_l==1
        title('Dropout')
    else
        title('College')
    end 
    hold on
    xticks(-50:10:50)
    set(gca,'FontName','Times New Roman','FontSize',FS-2);
    set(gcf,'color','w')
   end
    xlim([-50 50])



%% Descriptive wages across cohorts
clear all
close all
clc


% Read tuition data from Donovan and Herrington RED
cd('C:\Users\jbueren\Google Drive\endo_health\data\tuition_cost')
% Time college completion
Tcollege = 4;

% cohorts considered
start_year    = 1918;   % 1918
end_year      = 1990;   % 1990
year_grid     = start_year:end_year;
nt            = length(year_grid);

% Model preliminaries
R    = 1.04;             % gross return on svings

% CPI data
cpimat          = dlmread('cpi.txt');
cpi             = cpimat(find(cpimat(:,1)==start_year):find(cpimat(:,1)==end_year),2);  

% Public schools tuition and fees per FTE student
collegecost     = dlmread('ccosts_tf_public_fte.txt');
colcost_data    = (collegecost(find(collegecost(:,1) == start_year):find(collegecost(:,1) == (end_year+Tcollege-1)),2) .* cpimat(find(cpimat(:,1)==start_year):find(cpimat(:,1)==(end_year+Tcollege-1)),2))';

% PV of expected college cost
pv_colcost = zeros(1,nt);
for t = 1:nt
    pv_colcost(t) = sum( colcost_data(t:(t+3))./(R.^(0:3))   );
end


cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')

pattern = { '--'  '-'  '-.' ':' '-'};
colors =  {[0.0000, 0.4470, 0.7410];...   % Dark Blue
          [0.8500, 0.3250, 0.0980];...   % Dark Red
          [0.4660, 0.6740, 0.1880];...   % Dark Green
          [0.3, 0.1, 0.4]};  % Dark Purple
       

FS=10
figure(1)
set(1,'position',[150    150    450    250])
subplot(1,2,1)
for c_l=[ 2  4]
    data=importdata(strcat('mean_income_c_',int2str(c_l),'.txt'));
    plot(25:2:65,data(1:21,2)./2,'linewidth',1.5,'linestyle',pattern{c_l},'Color',colors{1})
    hold on
    plot(25:2:65,data(1:21,3)./2,'linewidth',1.5,'linestyle',pattern{c_l},'Color',colors{2})
    plot(25:2:65,data(1:21,4)./2,'linewidth',1.5,'linestyle',pattern{c_l},'Color',colors{3})
%     legend('Dropout 1930','College 1930','Dropout 1970','College 1970')
    ylim([0 90])
end
title('Mean Wage')
ylabel('$ (000s)')
xlabel('Age')
set(gca,'FontName','Times New Roman','FontSize',FS-2);
subplot(1,2,2)
plot(year_grid-18,pv_colcost./1000,'k')
xlabel('Birth year cohort')
xlim([start_year-18 end_year-18]);
ylim([0 12]);
title('Tuition costs')
set(gca,'FontName','Times New Roman','FontSize',FS-2);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\mean_wage_tuition','-depsc')


