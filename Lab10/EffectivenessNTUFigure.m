% Plot effectiveness versus NTU and Cr for crossflow heat exchanger

Cr=[0,.25,.5,.75];
NTU=linspace(0,1,100)';
e=zeros(100,5);

figure('color','w');
axes('position',[.15,.15,.75,.75],'nextplot','add','fontsize',16,...
    'xlim',[0,1],'ylim',[0,.7],'xminortick','on','yminortick','on',...
    'xgrid','on','ygrid','on','xminorgrid','on','yminorgrid','on',...
    'tickdir','out');
xlabel('NTU','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex','fontsize',18);
htitle=title('Effectiveness versus NTU for a crossflow heat exchanger',...
    'fontsize',16,'units','normalized','position',[.5,1.05]);

for (k=1:length(Cr))
    e(:,k)=(1-exp(-NTU*(1-Cr(k))))./(1-Cr(k)*exp(-NTU*(1-Cr(k))));
    hl(k)=plot(NTU,e(:,k),'k-','linewidth',1.5,'color',[.6,.6,.6]);
    ht(k)=text(NTU(end)*1.02,e(end,k),['$',num2str(Cr(k)),'$'],...
        'interpreter','latex','fontsize',16);
end

k=k+1;
e(:,k)=NTU./(1+NTU);
hl(k)=plot(NTU,e(:,k),'k-','linewidth',1.5,'color',[.6,.6,.6]);
ht(k)=text(NTU(end)*1.02,e(end,k),['$1$'],...
        'interpreter','latex','fontsize',16);
    
htCr=text(NTU(end)*1.04,e(end,1)*1.05,'$C_r$',...
    'interpreter','latex','fontsize',16,'horizontalalignment','center');