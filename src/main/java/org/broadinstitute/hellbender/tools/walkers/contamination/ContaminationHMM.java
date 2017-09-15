package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleBinaryOperator;

/**
 * PileupSummary: observed data
 * PileupSumary: the position class
 * Boolean: hidden state -- is is LoH or not
 */
public class ContaminationHMM implements HMM<PileupSummary, PileupSummary, Boolean> {
    // we consider allele fractions in a small range around 0.5 to be heterozygous.  Beyond that range is LoH.
    public static final double MIN_HET_AF = 0.3;
    public static final double MAX_HET_AF = 0.7;
    private static final double SNV_ERROR_RATE = 0.01;
    private static final double priorFractionOfGenomeWithLoH = 0.01;
    private static final double averageLoHSize = 5_000_000;
    private static final double averageNonLoHSize = 50_000_000;

    private final double contamination;

    public ContaminationHMM(final double contamination) {
        this.contamination = contamination;
    }

    public List<Boolean> hiddenStates() { return Arrays.asList(Boolean.TRUE, Boolean.FALSE); }

    public double logPriorProbability(final Boolean isLoH, final PileupSummary position) {
        return isLoH ? Math.log(priorFractionOfGenomeWithLoH) : Math.log(1 - priorFractionOfGenomeWithLoH);
    }

    public double logTransitionProbability(final Boolean fromLoH, final PileupSummary currentPosition, final Boolean toLoH, final PileupSummary nextPosition) {
        final double d = distance(currentPosition, nextPosition);
        final double probabilityOfSameState = Math.exp(-d / (fromLoH ? averageLoHSize : averageNonLoHSize));
        return Math.log(fromLoH == toLoH ? probabilityOfSameState : 1 - probabilityOfSameState);
    }

    public double logEmissionProbability(final PileupSummary ps, final Boolean isLoH, final PileupSummary position) {
        final double alleleFrequency = ps.getAlleleFrequency();
        final double homRefPrior = MathUtils.square(1 - alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homAltPrior = MathUtils.square(alleleFrequency);

        final int depth = ps.getTotalCount();
        final int altCount = ps.getAltCount();

        final DoubleBinaryOperator likelihood = (min, max) -> MathUtils.uniformBinomialProbability(depth, altCount, min, max);

        final double homRefLikelihood = likelihood.applyAsDouble(0, SNV_ERROR_RATE + contamination);
        final double homAltLikelihood = likelihood.applyAsDouble(1 - SNV_ERROR_RATE - contamination, 1);
        final double hetLikelihood = !isLoH ? likelihood.applyAsDouble(MIN_HET_AF, MAX_HET_AF) :
                0.5*likelihood.applyAsDouble(0, MIN_HET_AF) + 0.5*likelihood.applyAsDouble( MAX_HET_AF, 1);

        return Math.log(homRefPrior*homRefLikelihood + hetPrior*hetLikelihood + homAltPrior*homAltLikelihood);
    }

    private static double distance(final PileupSummary from, final PileupSummary to) {
        return from.getContig().equals(to.getContig()) ? Math.abs(to.getStart() - from.getStart()) : Double.POSITIVE_INFINITY;
    }
}
