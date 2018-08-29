function VBA_title(h,txt)
% VBA_title(h,txt)
% This function serves as a wrapper of the buil-in title() function that
% moreover ensures styling of the text
% IN: 
%   - h  : handle to the  plot to label  
%   - txt: text of the label  

t=title(h, txt);
set(t,'FontName'    ,'Arial'    );
set(t,'FontUnits'   ,'points'   );
set(t,'FontSize'    ,12         );
set(t,'FontWeight'  ,'bold'     );
set(t,'Color'       ,[0 0 0]    );
set(t,'Margin'      ,3          );