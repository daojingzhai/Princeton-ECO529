% This file makes 2D and 3D density diffusion movie for distributional
% impulse response. 
% Platform: MATLAB R2019a 
% Data require: solution_DIR.mat

clc;clear
load('solution_DIR.mat');
%% Create and open the video object
vidObj = VideoWriter('diffusion_animation_2D.avi');
open(vidObj);

% Open a figure and crate the axes
currfig = figure(8);
axes;

p1 = plot(Eta,pdf_stat,'r','LineWidth',2);
hold on

% Loop over the data to create the video
for i=1:size(sol1,1)
   % Plot the data
   h0 =plot(Eta,sol1(i,:),'Color',[1/3 1/3 1/3],'LineWidth',1);
   hold on
   legend([h0 p1],'Density Diffusion','Stationary Distribution','FontSize',14);
   str = ['t = ' num2str(T2(i),'%d')];
   h1 = text(.1,14,str,'Color','k','FontSize',16);
   xlabel('$\eta^e$','Interpreter','LaTex');
   ylabel('$f(\eta^e)$','Interpreter','LaTex');
   set(gca,'xlim',[0 1],'ylim',[0 18])
   % Get the current frame
   currFrame = getframe(currfig);
   % Write the current frame
   writeVideo(vidObj,currFrame);
   delete(h0)
   delete(h1)
end

p2 = plot(Eta,pdf_init,'b','LineWidth',2);
hold on

for i=1:size(sol2,1)
   % Plot the data
   h0 = plot(Eta,sol2(i,:),'Color',[1/3 1/3 1/3],'LineWidth',1);
   hold on
   legend([h0 p1 p2],'Density Diffusion','Stationary Distribution','Initial Distribution','FontSize',14);
   str = ['t = ' num2str(T2(i),'%d')];
   h1 = text(.1,14,str,'Color','k','FontSize',16);
   xlabel('$\eta^e$','Interpreter','LaTex');
   ylabel('$f(\eta^e)$','Interpreter','LaTex');
   set(gca,'xlim',[0 1],'ylim',[0 18])
   % Get the current frame
   currFrame = getframe(currfig);
   % Write the current frame
   writeVideo(vidObj,currFrame);
   delete(h0)
   delete(h1)
end

% slose (and save) the video object
close(vidObj);


%% 3D density diffusion movie
N = length(Eta);
Nt1 = length(T1);
Nt2 = length(T2);
Nt = Nt1+Nt2;

vidObj = VideoWriter('diffusion_animation_3D.avi');
open(vidObj);

% Open a figure and crate the axes
currfig = figure(9);
p0 = surf(Eta,fliplr(T0*20),sol0,'EdgeColor','none','FaceColor',[1/3,1/3,0],'FaceAlpha',0.3); hold on
p1 = plot3(Eta,ones(N)*(-T1(end)*20),pdf_stat,'LineWidth',2,'Color','r');
%  labels   
xlabel('$\eta^e$','Interpreter','LaTex');
ylabel('Time $t$','Interpreter','LaTex');
zlabel('Density $f(\eta^e,t)$','Interpreter','LaTex');
yticks([-50 -20 0 50 100 150 200])
yticklabels({'-5','-1','0','50','100','150','200'})
ylim([-50 200])
zlim([0 18])
set(gca,'Ydir','reverse')
text1 = text(0.8,100,15,{'Shock' 'Period'},'Color','magenta','FontSize',18);
% view(-30,21) 
% view(-3,23)
view(3,23)
% Loop over the data to create the video
for i=1:Nt-1
   if i < size(T1,2)
       h_temp_1 = surf(Eta,fliplr(-T1(end-i:end)*20),sol1(1:i+1,:),'EdgeColor','none','FaceColor','magenta','FaceAlpha',0.3);
       h_temp_2 = plot3(Eta,ones(N)*fliplr(-T1(end-i)*20), sol1(i+1,:),'LineWidth',2,'Color',[1/3 1/3 1/3]); 
   elseif i == size(T1,2)
       delete(text1)
       p2 = surf(Eta,fliplr(-T1*20),sol1,'EdgeColor','none','FaceColor','magenta','FaceAlpha',0.3); 
       p3 = plot3(Eta,ones(N)*T2(1),  pdf_init,'LineWidth',2,'Color','b'); 
       text2 = text(0.8,100,15,{'Transitional' 'Period'},'Color',[1/3,1/3,1/3],'FontSize',18); 
   elseif i > size(T1,2)
       h_temp_1 = surf(Eta,T2(1:i-size(T1,2)+1),sol2(1:i-size(T1,2)+1,:),'EdgeColor','none','FaceColor',[1/3,1/3,1/3],'FaceAlpha',0.3); 
       h_temp_2 = plot3(Eta,ones(N)*T2(i-size(T1,2)+1), sol2(i-size(T1,2)+1,:),'LineWidth',2,'Color',[1/3 1/3 1/3]); 
   end
   % Get the current frame
   currFrame = getframe(currfig);
   % Write the current frame
   writeVideo(vidObj,currFrame);
   delete(h_temp_1)
   delete(h_temp_2)
end
% Get the current frame
currFrame = getframe(currfig);
% Write the current frame
writeVideo(vidObj,currFrame);
% slose (and save) the video object
close(vidObj);