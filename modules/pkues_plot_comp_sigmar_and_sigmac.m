sigmac = (abs(Zp_norm).^2-abs(Zm_norm).^2)./(abs(Zp_norm).^2+abs(Zm_norm).^2);
sigmac(:,4,:,jpl) = (abs(Zp_norm(:,1,:,jpl)).^2+abs(Zp_norm(:,2,:,jpl)).^2 - ...
                     abs(Zm_norm(:,1,:,jpl)).^2-abs(Zm_norm(:,2,:,jpl)).^2)./...
                    (abs(Zp_norm(:,1,:,jpl)).^2+abs(Zp_norm(:,2,:,jpl)).^2 + ...
                     abs(Zm_norm(:,1,:,jpl)).^2+abs(Zm_norm(:,2,:,jpl)).^2);
sigmac(:,5,:,jpl) = (abs(Zp_norm(:,1,:,jpl)).^2+abs(Zp_norm(:,2,:,jpl)).^2+abs(Zp_norm(:,3,:,jpl)).^2 - ...
                     abs(Zm_norm(:,1,:,jpl)).^2-abs(Zm_norm(:,2,:,jpl)).^2-abs(Zm_norm(:,3,:,jpl)).^2)./...
                    (abs(Zp_norm(:,1,:,jpl)).^2+abs(Zp_norm(:,2,:,jpl)).^2+abs(Zp_norm(:,3,:,jpl)).^2 + ...
                     abs(Zm_norm(:,1,:,jpl)).^2+abs(Zm_norm(:,2,:,jpl)).^2+abs(Zm_norm(:,3,:,jpl)).^2);
for ixx=1:npa
    for s=1:S
        sigmar(ixx,1,s,jpl) = (abs(dVnorm(ixx,1,s,jpl)).^2-abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2)./(abs(dVnorm(ixx,1,s,jpl)).^2+abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2);
        sigmar(ixx,2,s,jpl) = (abs(dVnorm(ixx,2,s,jpl)).^2-abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2)./(abs(dVnorm(ixx,2,s,jpl)).^2+abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2);
        sigmar(ixx,3,s,jpl) = (abs(dVnorm(ixx,3,s,jpl)).^2-abs(Pola_norm(ixx,jpb,jpl,6)*sqrt(ns0(1)/ns0(s))).^2)./(abs(dVnorm(ixx,3,s,jpl)).^2+abs(Pola_norm(ixx,jpb,jpl,6)*sqrt(ns0(1)/ns0(s))).^2);
        sigmar(ixx,4,s,jpl) = (abs(dVnorm(ixx,1,s,jpl)).^2+abs(dVnorm(ixx,2,s,jpl)).^2 - ...
                               abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2-abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2)./...
                              (abs(dVnorm(ixx,1,s,jpl)).^2+abs(dVnorm(ixx,2,s,jpl)).^2 + ...
                               abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2+abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2);
        sigmar(ixx,5,s,jpl) = (abs(dVnorm(ixx,1,s,jpl)).^2+abs(dVnorm(ixx,2,s,jpl)).^2+abs(dVnorm(ixx,3,s,jpl)).^2 - ...
                               abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2-abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2-abs(Pola_norm(ixx,jpb,jpl,6)*sqrt(ns0(1)/ns0(s))).^2)./...
                              (abs(dVnorm(ixx,1,s,jpl)).^2+abs(dVnorm(ixx,2,s,jpl)).^2+abs(dVnorm(ixx,3,s,jpl)).^2 + ...
                               abs(Pola_norm(ixx,jpb,jpl,4)*sqrt(ns0(1)/ns0(s))).^2+abs(Pola_norm(ixx,jpb,jpl,5)*sqrt(ns0(1)/ns0(s))).^2+abs(Pola_norm(ixx,jpb,jpl,6)*sqrt(ns0(1)/ns0(s))).^2);
    end
end

if (S==3)
    h4 = figure('unit','normalized','Position',[0.01 0.1 0.8 0.8],...
      'DefaultAxesFontSize',15);
    subplot(351); hold on; box on;
    plot(pas,sigmar(:,1,1,1),'r',...
         pas,sigmac(:,1,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    legend('\sigma_r','\sigma_c');
    xlabel(strpa);ylabel('X');
    title('proton core')

    subplot(352); hold on; box on;
    plot(pas,sigmar(:,2,1,1),'r',...
         pas,sigmac(:,2,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Y');

    subplot(353); hold on; box on;
    plot(pas,sigmar(:,3,1,1),'r',...
         pas,sigmac(:,3,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Z');

    subplot(354); hold on; box on;
    plot(pas,sigmar(:,4,1,1),'r',...
         pas,sigmac(:,4,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('perp');

    subplot(355); hold on; box on;
    plot(pas,sigmar(:,5,1,1),'r',...
         pas,sigmac(:,5,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('trace');


    subplot(356); hold on; box on;
    plot(pas,sigmar(:,1,2,1),'r',...
         pas,sigmac(:,1,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    legend('\sigma_r','\sigma_c');
    xlabel(strpa);ylabel('X');
    title('proton beam')

    subplot(357); hold on; box on;
    plot(pas,sigmar(:,2,2,1),'r',...
         pas,sigmac(:,2,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Y');

    subplot(358); hold on; box on;
    plot(pas,sigmar(:,3,2,1),'r',...
         pas,sigmac(:,3,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Z');

    subplot(359); hold on; box on;
    plot(pas,sigmar(:,4,2,1),'r',...
         pas,sigmac(:,4,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('perp');

    subplot(3,5,10); hold on; box on;
    plot(pas,sigmar(:,5,2,1),'r',...
         pas,sigmac(:,5,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('trace');

    subplot(3,5,11); hold on; box on;
    plot(pas,sigmar(:,1,3,1),'r',...
         pas,sigmac(:,1,3,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    legend('\sigma_r','\sigma_c');
    xlabel(strpa);ylabel('X');
    title('electron')

    subplot(3,5,12); hold on; box on;
    plot(pas,sigmar(:,2,3,1),'r',...
         pas,sigmac(:,2,3,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Y');

    subplot(3,5,13); hold on; box on;
    plot(pas,sigmar(:,3,3,1),'r',...
         pas,sigmac(:,3,3,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Z');

    subplot(3,5,14); hold on; box on;
    plot(pas,sigmar(:,4,3,1),'r',...
         pas,sigmac(:,4,3,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('perp');

    subplot(3,5,15); hold on; box on;
    plot(pas,sigmar(:,5,3,1),'r',...
         pas,sigmac(:,5,3,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('trace');
end

if (S==2)
    h4 = figure('unit','normalized','Position',[0.01 0.1 0.8 0.6],...
      'DefaultAxesFontSize',15);
    subplot(251); hold on; box on;
    plot(pas,sigmar(:,1,1,1),'r',...
         pas,sigmac(:,1,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    legend('\sigma_r','\sigma_c');
    xlabel(strpa);ylabel('X');
    title('proton')

    subplot(252); hold on; box on;
    plot(pas,sigmar(:,2,1,1),'r',...
         pas,sigmac(:,2,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Y');

    subplot(253); hold on; box on;
    plot(pas,sigmar(:,3,1,1),'r',...
         pas,sigmac(:,3,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Z');

    subplot(254); hold on; box on;
    plot(pas,sigmar(:,4,1,1),'r',...
         pas,sigmac(:,4,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('perp');

    subplot(255); hold on; box on;
    plot(pas,sigmar(:,5,1,1),'r',...
         pas,sigmac(:,5,1,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('trace');


    subplot(256); hold on; box on;
    plot(pas,sigmar(:,1,2,1),'r',...
         pas,sigmac(:,1,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    legend('\sigma_r','\sigma_c');
    xlabel(strpa);ylabel('X');
    title('electron')

    subplot(257); hold on; box on;
    plot(pas,sigmar(:,2,2,1),'r',...
         pas,sigmac(:,2,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Y');

    subplot(258); hold on; box on;
    plot(pas,sigmar(:,3,2,1),'r',...
         pas,sigmac(:,3,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('Z');

    subplot(259); hold on; box on;
    plot(pas,sigmar(:,4,2,1),'r',...
         pas,sigmac(:,4,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('perp');

    subplot(2,5,10); hold on; box on;
    plot(pas,sigmar(:,5,2,1),'r',...
         pas,sigmac(:,5,2,1),'b','linewidth',2);
    xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
    ylim([-1,1]);
    xlabel(strpa);ylabel('trace');
end

print(h4,'-dpng',[savepath,'fig_pdrk_',figstr,'_sigmar_sigmac.png']);
savefig([savepath,'fig_pdrk_',figstr,'_sigmar_sigmac.fig']);