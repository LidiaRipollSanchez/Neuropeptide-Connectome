function [fractionAC_allNeurons,fractionAC_Pharynx, fractionAC_Sensory, fractionAC_Inter, fractionAC_Motor] = ComputeFractionAutocrineConnections(CoExpressionMatrix, anatomical_class)
% Jan Watteyne
% function to calculate the fraction of neurons showing putative autocrine
% connections (~co-expression of both NPP and GPCR). Both for all neurons
% as specified for pharyngeal - sensory neurons - interneurons and
% motoneurons

neuronWithOrWithoutCoExpression = sum(CoExpressionMatrix, 1);
numberNeuronsWithCoExpression = length(find(neuronWithOrWithoutCoExpression ~= 0));
fractionAC_allNeurons = numberNeuronsWithCoExpression / (size(CoExpressionMatrix, 2)-1);

indexPharynx = find(contains(anatomy.FinalClassification,'Pharynx'));
neuronWithOrWithoutCoExpression_Pharynx = neuronWithOrWithoutCoExpression(indexPharynx);
fractionAC_Pharynx = length(find(neuronWithOrWithoutCoExpression_Pharynx ~= 0)) / (length(neuronWithOrWithoutCoExpression_Pharynx));

indexSensory = find(contains(anatomy.FinalClassification,'sensory neuron'));
neuronWithOrWithoutCoExpression_Sensory = neuronWithOrWithoutCoExpression(indexSensory);
fractionAC_Sensory = length(find(neuronWithOrWithoutCoExpression_Sensory ~= 0)) / (length(neuronWithOrWithoutCoExpression_Sensory));

indexInter = find(contains(anatomy.FinalClassification,'interneuron'));
neuronWithOrWithoutCoExpression_Inter = neuronWithOrWithoutCoExpression(indexInter);
fractionAC_Inter = length(find(neuronWithOrWithoutCoExpression_Inter ~= 0)) / (length(neuronWithOrWithoutCoExpression_Inter));

indexMotor = find(contains(anatomy.FinalClassification,'motor neuron'));
neuronWithOrWithoutCoExpression_Motor = neuronWithOrWithoutCoExpression(indexMotor);
fractionAC_Motor = length(find(neuronWithOrWithoutCoExpression_Motor ~= 0)) / (length(neuronWithOrWithoutCoExpression_Motor));
end

