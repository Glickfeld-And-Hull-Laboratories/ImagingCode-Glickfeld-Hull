function fix_axes(fig_name,font_size,x_label,y_label)

figure(fig_name)
hline = findobj(gcf,'type','line');
set(gca,'TickDir','out','FontSize',font_size,'Color','none')
%set(hline,'LineWidth',1.2);
box(gca,'off');
xlabel(x_label)
ylabel(y_label)
%axis square
end