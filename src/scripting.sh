#for folder in fourKnotLsAnnealingSpan/ fourKnotLsAnnealingSpanRefin/ fourKnotLsAnnealingTangent/ fourKnotLsAnnealingTangentRefin/ fourKnotLsGeneticSpan/ fourKnotLsGeneticTangent/ needleLsAnnealingSpan/ needleLsAnnealingTangent/ needleLsGeneticSpan/ needleLsGeneticTangent/
for folder in fourKnotLsAnnealingTangent/ fourKnotLsAnnealingTangentRefin/ fourKnotLsGeneticTangent/ needleLsAnnealingTangent/  needleLsGeneticTangent 
do
    julia lsAnalysisDistributed.jl --path ../outputs/$folder --burnout 0.95 --name  max_span_minimum --minFunc squaredMaxSpan  --processes 40 --minimum
    #cat ../outputs/${folder}/metaParams.txt
    echo ""
	echo "done"
    echo ""
done
