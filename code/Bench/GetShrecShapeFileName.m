function shapeFileName = GetShrecShapeFileName(number,class,strength)
shapeFileName = [sprintf('%04d',number) '.' class '.' int2str(strength)];