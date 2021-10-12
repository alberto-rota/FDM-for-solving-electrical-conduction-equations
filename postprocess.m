function postprocess(fname)
%
% Plot results from the simultion performed
% with fd1D_FK.m
%
% fname: filename to postprocess (generated by fd1D_FK), e.g., output.mat
%
close all;
load(fname);
%
% Node number for output
%
nnd = length(x);
dx = x(2)-x(1);
DN = ceil(0.5/dx);
nodeOut=[ceil(nnd/2)-DN ceil(nnd/2) ceil(nnd/2)+DN];
nnout=length(nodeOut);
CV = CV_computation(x, Vsol);
save(fname,'CV','-append');
% plot_AP(x, nodeOut, Vsol(:,1), Vsol(:,nodeOut+1));
% plot_propagation(x,Vsol,Jion);
end

function CV = CV_computation(node,Vm)
%
% CV is calculated based in the average of 4 segments within the cable
%
nnd = length(node);
dL = 0.5;  % length to calculate velocity
dx = node(2)-node(1);
spc = fix(dL/dx); % number of nodes present in the interval

midNd = fix(nnd/2);
Nd = [midNd-2*spc midNd-spc midNd midNd+spc midNd+2*spc];

[Y,I]=max(diff(Vm(:,Nd(1)+1))./diff(Vm(:,1)));
t1=Vm(I,1);
if(Y<0.2)
    fprintf('Conduction velocity cannot be computed reliably\n');
    fprintf('Increase simulation time ...\n');
    return;
end
cont = 1;
CV=0.0;
for i=1:4
  if(max(Vm(:,Nd(i+1)+1))<1.0)
    fprintf('The front is not reaching the end of the cable...\n');
    fprintf('Increase simulation time to improve CV estimation\n');
    fprintf('Conduction velocity computed with %d measurements\n',cont-1);
    break;
  end
  [Y,I]=max(diff(Vm(:,Nd(i+1)+1))./diff(Vm(:,1)));
  t2=Vm(I,1);
  if((Y<0.2)&&(cont<2))
    fprintf('Conduction velocity cannot be computed reliably\n');
    fprintf('Increase simulation time ...\n');
    return;
  else
      dx = node(Nd(i+1))-node(Nd(i));
      CV1 = 1000*dx/(t2-t1);
      CV = CV + CV1;
      cont = cont+1;
      t1=t2;
  end
end
CV = CV/(cont-1);
fprintf('Conduction velocity (cm/s): %6.2f\n',CV);
end

function plot_AP(node, nodeOut, t, AP)
%
nnd = length(node);
%
% Drawing Geometry
%
figure1 = figure;

% Create subplot
subplot1 = subplot(2,2,1,'Parent',figure1,'YColor',[1 1 1],'FontSize',10,...
    'YTick',[]);
ylim(subplot1,[-0.5 1]);
hold(subplot1,'on');
% Create plot
plot(node,zeros(nnd,1),'DisplayName','Geometry','LineWidth',2,'Color',[0 0 0]);
% Create plot
plot(node(nodeOut(1)),0,'DisplayName','Cell 1','MarkerFaceColor',[0 0.498039215803146 0],...
    'MarkerEdgeColor',[0 0.498039215803146 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0.498039215803146 0]);
% Create plot
plot(node(nodeOut(2)),0,'DisplayName','Cell 2','MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
% Create plot
plot(node(nodeOut(3)),0,'DisplayName','Cell 3','MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
% Create xlabel
xlabel('Position (cm)');
% Create title
title('Geometry');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,'FontSize',10,'Location','northwest');
%
% Drawing Action Potentials
%
% Create subplot
subplot2 = subplot(2,2,2,'Parent',figure1,'FontSize',10,...
    'YTickLabel',{'0.0','0.5','1.0','1.5'},...
    'YTick',[0 0.5 1 1.5]);
ylim(subplot2,[-0.1 1.6]);
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'on');
% Create plot
plot(t,AP(:,1),'Parent',subplot2,'DisplayName','Cell 1',...
    'LineWidth',2,'Color',[0 0.498039215803146 0]);
xlabel('Time (ms)');
ylabel('Normalized Potential (-)');
% Create title
title('Cell 1');
%
% Create subplot
subplot3 = subplot(2,2,3,'Parent',figure1,'FontSize',10,...
    'YTickLabel',{'0.0','0.5','1.0','1.5'},...
    'YTick',[0 0.5 1 1.5]);
ylim(subplot3,[-0.1 1.6]);
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'on');
% Create plot
plot(t,AP(:,2),'Parent',subplot3,'DisplayName','Cell 1',...
    'LineWidth',2,'Color',[1 0 0]);
xlabel('Time (ms)');
ylabel('Normalized Potential (-)');
% Create title
title('Cell 2');
%
% Create subplot
subplot4 = subplot(2,2,4,'Parent',figure1,'FontSize',10,...
    'YTickLabel',{'0.0','0.5','1.0','1.5'},...
    'YTick',[0 0.5 1 1.5]);
ylim(subplot4,[-0.1 1.6]);
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'on');
% Create plot
plot(t,AP(:,3),'Parent',subplot4,'DisplayName','Cell 1',...
    'LineWidth',2,'Color',[0 0 1]);
xlabel('Time (ms)');
ylabel('Normalized Potential (-)');
% Create title
title('Cell 3');
end

function plot_propagation(node,Vm,Jion)

nnd = length(node);
nstep = length(Vm(:,1));

figure2=figure;
subplot1 = subplot(1,2,1,'Parent',figure2,'FontSize',16,...
    'YTickLabel',{'0.0','0.5','1.0','1.5'},'YTick',[0 0.5 1 1.5],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 100 200 300 400]);
   xlabel('Position (cm)');
   ylabel('Normalized Potential (-)');
   ylim(subplot1,[-0.1 1.6]);
   drawnow;
subplot2 = subplot(1,2,2,'Parent',figure2,'FontSize',16,...
    'YTickLabel',{'-2.5','-2.0','-1.5','-1.0','-0.5','0.0'},...
    'YTick',[-2.5,-2.0 -1.5 -1.0 -0.5 0.0],...
    'XTickLabel',{'0','100','200','300','400'},...
    'XTick',[0 100 200 300 400]);
   xlabel('Position (cm)');
   ylabel('Ionic Current (-)');
   ylim(subplot2,[-2.6 0.1]);
   drawnow;
flag = 1;
while(flag)
 for i=1:nstep
   plot(node,Vm(i,2:end),'Parent',subplot1,'LineWidth',2,'Color',[0 0 0]);
   ylim(subplot1,[-0.1 1.6]);
   drawnow;
   plot(node,Jion(i,2:end),'Parent',subplot2,'LineWidth',2,'Color',[0 0 0]);
   ylim(subplot2,[-2.6 0.1]);
   drawnow;
 end
 flag = input('Repeat (0/1): ');
end
end