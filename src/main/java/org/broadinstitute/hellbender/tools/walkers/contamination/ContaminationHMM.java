package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.hmm.HMM;

import java.util.Arrays;
import java.util.List;

/**
 * PileupSummary: observed data
 * PileupSumary: the position class
 * Boolean: hidden state -- is is LoH or not
 */
public class ContaminationHMM implements HMM<PileupSummary, PileupSummary, Boolean> {
    // we consider allele fractions in a small range around 0.5 to be heterozygous.  Beyond that range is LoH.
    public static final double MIN_HET_AF = 0.4;
    public static final double MAX_HET_AF = 0.6;
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
        final double hetMin = MIN_HET_AF - contamination / 2;
        final double hetMax = MAX_HET_AF + contamination / 2;

        final double homRefLikelihood = MathUtils.uniformBinomialProbability(depth, altCount, 0, contamination);
        final double homAltLikelihood = MathUtils.uniformBinomialProbability(depth, altCount, 1 - contamination, 1);
        final double hetLikelihood = isLoH ? MathUtils.uniformBinomialProbability(depth, altCount, contamination, hetMin) +
                MathUtils.uniformBinomialProbability(depth, altCount, hetMax, 1 - contamination) :
                MathUtils.uniformBinomialProbability(depth, altCount, hetMin, hetMax);

        return Math.log(homRefPrior*homRefLikelihood + hetPrior*hetLikelihood + homAltPrior*homAltLikelihood);

    }

    private static double distance(final PileupSummary from, final PileupSummary to) {
        return from.getContig().equals(to.getContig()) ? Math.abs(to.getStart() - from.getStart()) : Double.POSITIVE_INFINITY;
    }
}
