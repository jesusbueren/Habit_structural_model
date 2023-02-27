clear all
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')
fileID = fopen('asset_distribution_model.txt','r');
A=fscanf(fileID,'%f');
fclose(fileID);
A = standardizeMissing(A,-9)
T=36
type_y=3
type_e=3
A=reshape(A,7,T,type_y,type_e)

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5]
FS=10
figure(1)
set(1,'position',[150    150    750    750])
for p_l=1:3
for e_l=1:type_e
    subplot(3,3,(p_l-1)*3+e_l)
    for y_l=1:type_y
        if e_l==1
            h(y_l)=plot(26:2:80,A(3+p_l,1:28,y_l,e_l),'Color',colors{y_l},'linewidth',1.5,'linestyle',pattern{y_l})
        else
            plot(26:2:80,A(3+p_l,1:28,y_l,e_l),'Color',colors{y_l},'linewidth',1.5,'linestyle',pattern{y_l})
        end 
        hold on
    end 
    if e_l==1
        title('dropout')
%         ylim([0 150])
    elseif e_l==2
        title('highschool')
%         ylim([0 600])
    else
        title('college')
%         ylim([0 1000])
    end
    if e_l==1
        if p_l==1
            ylabel('P25')
        elseif p_l==2
            ylabel('P50')
        elseif p_l==3
            ylabel('P75')
        end 
    end
ylim([0 2000])
    set(gca,'FontName','Times New Roman','FontSize',FS);
end
end
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')

%% Costs
clear all
close all
clc
cd('C:\Users\jbueren\Google Drive\endo_health\structural_model')
data=importdata('costs_e_y.txt')

colors = {  [0.4660    0.6740    0.1880]+0.1   [0.9290    0.6940    0.1250]+0.1    [0.8500    0.3250    0.0980]+0.1 };
colors2 ={  [0.4660    0.6740    0.1880]-0.2   [0.9290    0.6940    0.1250]-0.2    [0.8500    0.3250    0.0980]-0.2   };
colors3 =   [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; 0.8500    0.3250    0.0980] ;
FS=10

data_mat=data.data
data_mat=reshape(data_mat(:,2:5),3,3,4)
costs=data_mat(:,:,1)
data_mat(:,:,1)=data_mat(:,:,1)-mean(costs(:))

for j=1:3
for e_l=1:3
figure(j)
set(j,'position',[150    150    750    325])
subplot(1,3,e_l)
X=categorical({'Protective','Detrimental','Harmful'})
X = reordercats(X,{'Protective','Detrimental','Harmful'});
if j<3
    b=bar(X,data_mat(:,e_l,j),'FaceColor','flat')
    for c_l=1:3
        b.CData(c_l,:) = colors{c_l}
    end
else
    b=bar(X,squeeze(data_mat(:,e_l,j:j+1))./sum(squeeze(data_mat(:,e_l,j:j+1))),'FaceColor','flat')
    for c_l=1:3
       b(1).CData(c_l,:) = colors{c_l}
       b(2).CData(c_l,:) = colors2{c_l}
    end
end

if e_l==1
    title('dropouts')
elseif e_l==2
    title('Highschool')
else
    title('college')
end
set(gca,'FontName','Times New Roman','FontSize',FS);
if j==2
    ylim([0,4.0])
elseif j==1
    ylim([-2,2])
elseif j==3
    ylim([0,1])
end
set(gcf,'color','w')
end
end

sz = 50;
x=data_mat(:,:,2)
y=log(data_mat(:,:,3))
figure(4)
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
s=scatter(x(:,1),y(:,1),sz,colors3,'filled','d')
s=scatter(x(:,2),y(:,2),sz+10,colors3,'filled','s')
s=scatter(x(:,3),y(:,3),sz,colors3,'filled','o')
ylim([-5 0])
xlim([1.2 4.0])
xlabel('Value function at age 25')
ylabel('log share')
theString = sprintf('y = %.3f x  %.3f ', Slope, Intercept);
text(3.2, -2.2, theString,'FontName','Times New Roman', 'FontSize', FS)
theString = sprintf('R^2= %.3f',Rsq);
text(3.2, -2.5, theString,'FontName','Times New Roman', 'FontSize', FS)
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')




