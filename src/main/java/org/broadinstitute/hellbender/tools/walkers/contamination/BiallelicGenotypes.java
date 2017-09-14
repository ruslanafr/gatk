package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.EnumMap;

public enum BiallelicGenotypes {
    HOM_REF, HET, HOM_ALT;

    public static class Posterior extends EnumMap<BiallelicGenotypes, Double> {
        public Posterior(final double homRefProb, final double hetProb, final double homAltProb) {
            super(BiallelicGenotypes.class);
            put(HOM_REF, homRefProb);
            put(HET, hetProb);
            put(HOM_ALT, homAltProb);
        }
    }

    public static Posterior getPosteriors(final PileupSummary ps, final double contamination) {
        final double alleleFrequency = ps.getAlleleFrequency();
        final double homRefPrior = MathUtils.square(1 - alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homAltPrior = MathUtils.square(alleleFrequency);

        final int totalCount = ps.getTotalCount();
        final int altCount = ps.getAltCount();

        final double homRefLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, 0, contamination);
        final double homAltLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, 1 - contamination, 1);
        final double hetLikelihood = MathUtils.uniformBinomialProbability(totalCount, altCount, ContaminationHMM.MIN_HET_AF - contamination / 2, ContaminationHMM.MAX_HET_AF + contamination / 2);

        final double[] unnormalized = new double[BiallelicGenotypes.values().length];
        unnormalized[BiallelicGenotypes.HOM_REF.ordinal()] = homRefLikelihood * homRefPrior;
        unnormalized[BiallelicGenotypes.HET.ordinal()] = hetLikelihood * hetPrior;
        unnormalized[BiallelicGenotypes.HOM_ALT.ordinal()] = homAltLikelihood * homAltPrior;
        final double[] normalized = MathUtils.normalizeFromRealSpace(unnormalized, true);

        return new Posterior(normalized[0], normalized[1], normalized[2]);
    }
}
