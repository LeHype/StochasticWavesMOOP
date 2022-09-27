
f = gcf;
ax = findobj(f, 'type', 'axes');
set(f,'color','w');
for i = 1:length(ax)
f = ax(i);
Linewidth = 6;
lines = findobj(f,'Type','Line','-or','Type','ErrorBar','-or','Type','Scatter');
for i = 1:numel(lines)  
  lines(i).LineWidth = Linewidth;
  if strcmp(lines(i).Type,  'scatter')
      lines(i).SizeData = 90;
  end
end

f.FontSize =40
try
f.Legend.FontSize = 40;
end


end
%%
f.CurrentAxes.XLabel.String='Number Components PLSR';
f.CurrentAxes.YLabel.String ='Root Mean squared Deviation';

%%
f = gcf;
Linewidth = 3;
for i = 1:numel(lines)
  lines(i).LineWidth = Linewidth;

end
f.CurrentAxes.FontSize =30
f.CurrentAxes.Legend.FontSize = 60;
set(f,'color','w');
lines(4).Color='#0072BD'
lines(4).LineStyle = ':'
lines(4).LineWidth= 8

lines(3).Color='#0072BD'
lines(3).LineStyle = ':'
lines(3).LineWidth= 8

% 
lines(2).Color = '#00cc00'
lines(2).LineStyle = '-'
lines(1).Color = 	'#cc0000'
lines(1).LineStyle = '-'


%%
f.CurrentAxes.XLabel.String='Cycles';
f.CurrentAxes.YLabel.String ='Target';
    

