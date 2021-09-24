%-------------------------------------------------------------------------
% ME EN 4650
%
% Drag and lift coefficients of the NACA 0012 airfoil
%
% Data from R.E. Sheldahl and P. C. Klimas, "Aerodynamic characteristics 
% of seven symmetrical airfoil sections through 180-degree angle of attack 
% for use in the aerodynamic analysis of vertical axis wind turbines", 
% Sandia National Laboratories Report, SAND80-2114, 1981
%
% M Metzger
% 3-2019
%-------------------------------------------------------------------------

% filename containing drag and lift coefficient as a function of angle 
% of attack and chord Reynolds number
FileName = 'NACA0012_SheldahlKlimas1981.csv';

% angle of attack values investigated
AoA = [0:25];  %deg

% Chord Reynolds numbers investigated
Rec = [1e4,2e4,4e4,8e4,1.6e5,3.6e5,7e5,1e6,2e6,5e6];
NumRec=length(Rec);

% read in data from file: delimiter=','; no. of header lines=5
Dat=importdata(FileName,',',5);

% extract the lift coefficient data
Cl=Dat.data(:,1:2:end);
Cd=Dat.data(:,2:2:end);

% create fiures for lift and drag
hf1=figure('color','w');
axes('position',[.12,.15,.7,.775],'nextplot','add','fontsize',16,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on',...
    'xminorgrid','on','yminorgrid','on','xlim',[0,20],'tickdir','out');

xlabel('$\alpha$ (deg)','interpreter','latex');
ylabel('$C_L$','interpreter','latex');

hf2=figure('color','w');
axes('position',[.12,.15,.7,.775],'nextplot','add','fontsize',16,...
    'xgrid','on','ygrid','on','xminortick','on','yminortick','on',...
    'xminorgrid','on','yminorgrid','on','xlim',[0,20],'tickdir','out');

xlabel('$\alpha$ (deg)','interpreter','latex');
ylabel('$C_D$','interpreter','latex');

% line colors
LC=fliplr(linspace(0,.9,NumRec));
LegStrCl=cell(NumRec);
LegStrCd=cell(NumRec);

% loop through the different chord Reynolds numbers
for (k=1:NumRec)
    % legend string
    LegStr{k}=num2str(Rec(k),'%3.1e');
    
    % plot lift coefficient
    figure(hf1);
    hlift(k)=plot(AoA,Cl(:,k),'-','color',ones(1,3)*LC(k),'linewidth',1.5);
    
    % plot drag coefficient
    figure(hf2);
    hdrag(k)=plot(AoA,Cd(:,k),'-','color',ones(1,3)*LC(k),'linewidth',1.5);
end

% add legends
figure(hf1);
hg1=legend(LegStr,'box','off','fontsize',12,'position',[.835,.277,.152,.365]);
annotation('textbox','string','$Re_c$','interpreter','latex',...
    'position',[.915,.625,.08,.065],'fontsize',14,'linestyle','none');

figure(hf2);
hg2=legend(LegStr,'box','off','fontsize',12,'position',[.835,.277,.152,.365]);
annotation('textbox','string','$Re_c$','interpreter','latex',...
    'position',[.915,.625,.08,.065],'fontsize',14,'linestyle','none');