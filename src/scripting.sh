#for folder in fourKnotLsAnnealingSpan/ fourKnotLsAnnealingSpanRefin/ fourKnotLsAnnealingTangent/ fourKnotLsAnnealingTangentRefin/ fourKnotLsGeneticSpan/ fourKnotLsGeneticTangent/ needleLsAnnealingSpan/ needleLsAnnealingTangent/ needleLsGeneticSpan/ needleLsGeneticTangent/
for folder in fourKnotLsAnnealingTangent/ fourKnotLsAnnealingTangentRefin/ fourKnotLsGeneticTangent/ needleLsAnnealingTangent/  needleLsGeneticTangent 
do
    julia lsAnalysisDistributed.jl --path ../../buenasSimulaciones/$folder --burnout 0.95 --name  max_span_mean --minFunc squaredMaxSpan  --processes 4
    #cat ../outputs/${folder}/metaParams.txt
    echo ""
	echo "done"
    echo ""
done
