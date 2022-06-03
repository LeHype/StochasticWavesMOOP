function poIdc = ParetoFilter(costValues)

eliminated = [];

for iCostPair = 1:size(costValues,1)-1
   
    compareWith = setdiffInt((iCostPair+1):size(costValues, 1), eliminated);
    
    dominationTest = costValues(iCostPair, :) - costValues(compareWith, :);
    
    dominated = any(all(sign(dominationTest) > -1,2), 1);
    if dominated
        eliminated = [eliminated, iCostPair];
    end
   
    dominates = all(sign(dominationTest) < 1, 2);
    if any(dominates)
        eliminated = [eliminated, compareWith(dominates)];
    end
end

poIdc = setdiffInt(1:size(costValues,1), eliminated);

end

