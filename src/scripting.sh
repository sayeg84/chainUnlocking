for folder in fourKnotLsAnnealingSpan/ fourKnotLsAnnealingSpanRefin/ fourKnotLsAnnealingTangent/ fourKnotLsAnnealingTangentRefin/ fourKnotLsGeneticSpan/ fourKnotLsGeneticTangent/ needleLsAnnealingSpan/ needleLsAnnealingTangent/ needleLsGeneticSpan/ needleLsGeneticTangent/ 
do
    #julia lsAnalysis.jl --path ../outouts/$folder --burnout 0.95 
    ls ../outouts/$folder/metaParams.txt
    echo "done"
done