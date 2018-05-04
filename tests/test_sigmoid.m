x = -2 : 2;

opt = struct('slope',3); 
err.slope = norm(x - VBA_sigmoid(VBA_sigmoid(x,opt),opt,'inverse',true));

opt = struct('center',3); 
err.center = norm(x - VBA_sigmoid(VBA_sigmoid(x,opt),opt,'inverse',true));

opt = struct('offset',3); 
err.offset = norm(x - VBA_sigmoid(VBA_sigmoid(x,opt),opt,'inverse',true));


opt = struct('scale',3); 
err.scale = norm(x - VBA_sigmoid(VBA_sigmoid(x,opt),opt,'inverse',true));

opt = struct('lapseRate',.3); 
err.lapse = norm(x - VBA_sigmoid(VBA_sigmoid(x,opt),opt,'inverse',true));



err



[~, d1]=VBA_sigmoid(-2:2,'slope',3); [~, d2]=sigm(-2:2,[],log(3)); norm(d1-d2)
[~, d1]=VBA_sigmoid(-2:2,'center',3); [~, d2]=sigm(-2:2,[],[0 3]'); norm(d1-d2)
[~, d1]=VBA_sigmoid(-2:2,'offset',3); [~, d2]=sigm(-2:2,struct('S0',3)); norm(d1-d2)
[~, d1]=VBA_sigmoid(-2:2,'scale',3); [~, d2]=sigm(-2:2,struct('G0',3)); norm(d1-d2)



