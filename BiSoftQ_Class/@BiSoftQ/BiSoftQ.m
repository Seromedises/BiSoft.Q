classdef BiSoftQ < handle

    properties
        
        PleatedMembrane PleatedMembrane
        Bellows Bellows
        type char
        
        p_gauge_MPa
        x_vec
        F_pull_vec
        F_push_vec
        F_el_vec
        F_el_membr_vec
        F_el_bell_vec
        VolumeChamber1_vec 
        VolumeChamber2_vec

        NPoints

        % effective strokes
        xtot_eff
        xmax_eff
        xmin_eff

        % efficiency
        deltaM
        Win
        Wout
        
    end
   
    methods (Access=public)

        function obj = BiSoftQ(PleatedMembraneObj,BellowsObj,type,p_gauge_MPa)

            arguments
               PleatedMembraneObj PleatedMembrane
               BellowsObj Bellows
               type (1,1) char {mustBeMember(type,{'A','B'})}
               p_gauge_MPa (1,1) double {mustBePositive}
            end

            % Constructor
            obj.PleatedMembrane = PleatedMembraneObj;
            obj.Bellows = BellowsObj;
            obj.type = type;
            obj.p_gauge_MPa = p_gauge_MPa;

        end
        
        AnimateActuator(obj,options)
        ForceCharacteristic(obj,options)
        PlotForceCharacteristic(obj)
        PlotNominalGeometry(obj)
        PlotVolumeChambers(obj)
        SummaryDesignParameters(obj,saveSummary,SummaryName)
        VolumeChambers(obj,options)
        Efficiency(obj,options)

    end


end