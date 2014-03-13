function grapher_staticDcmDisplay(nodes,node_colors,connect_colors,fig)


if nargin < 4
fig=figure('OuterPosition',[500 300 350 350],...
            'name','DCM Dsiplay',...
            'MenuBar','none','NumberTitle','off');  
end
hold on


if isempty(node_colors)
    node_colors = -10*eps*ones(1,nodes(1).n);
end

state = [];
for i=1:nodes(1).n
    state(i) = grapher_drawNode(nodes(i),node_colors(i));  
end

connections = [];
for i=1:nodes(1).n
    for j=1:nodes(1).n
        if i~=j
        if ~isnan(connect_colors(i,j)) % && connect_colors(i,j)~=0
        connections{i,j} = grapher_drawConnect(nodes(i),nodes(j),connect_colors(i,j));
        else
        connections{i,j} = NaN;
        end
        end
    end
end



% colormap(flipud(cbrewer('div','RdBu',100)));
colormap(my_dcm_colormap);
xlim([-200 200])
ylim([-200 200])
caxis([-1 1]);
% box on
axis off
set(gcf, 'color', [1 1 1])

drawnow



end

function c=my_dcm_colormap()

REDS = [...
	156 031 032;
	217 033 038;
	232 039 037;
	238 083 035;
	246 140 050;
	252 195 057;
	252 239 075];

BLUES = [...
	063 088 158;
	063 101 162;
	072 141 193;
	091 187 225;
	136 204 225;
	156 218 230;
	173 223 230;
	208	236	238];	

BLUES = cbrewer('seq','Blues',50);
REDS = cbrewer('seq','Reds',50);

c = [flipud(BLUES) ; REDS ];
end

