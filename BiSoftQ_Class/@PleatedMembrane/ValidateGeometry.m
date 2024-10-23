function [] = ValidateGeometry(obj,GeomDefWarnMsg)

if strcmp(GeomDefWarnMsg,"")
    GeomDefWarnMsg = lastwarn;
end

if ~strcmp(GeomDefWarnMsg,"")
    obj.ValidGeometry = false;
    if ~obj.ParametricStudyBool
       error(GeomDefWarnMsg) 
    end
else
    obj.ValidGeometry = true;
end

end