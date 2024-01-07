function txt = myupdatefcn(obj, event_obj)
pos = event_obj.Position;
txt = {sprintf('X: %.16g', pos(1)), sprintf('Y: %.16g', pos(2))};