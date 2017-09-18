package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.stream.IntStream;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static org.broadinstitute.hellbender.utils.MathUtils.log10Factorial;
import static org.broadinstitute.hellbender.utils.MathUtils.log10ToLog;

/**
 * Contains likelihood methods for the allele-fraction model.
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionLikelihoods {
//    private static final double LN_2 = Math.log(2.);
//
//    private AlleleFractionLikelihoods() {}
//
//    public static double hetLogLikelihood(final AlleleFractionGlobalParameters parameters,
//                                          final double minorFraction,
//                                          final AllelicCount count,
//                                          final AlleleFractionIndicator indicator) {
//        return hetLogLikelihood(parameters, minorFraction, count, indicator, AllelicPanelOfNormals.EMPTY_PON);
//    }
//
//    /**
//     * Computes {@link AlleleFractionLikelihoods#hetLogLikelihood} using allelic-bias priors derived from an
//     * {@link AllelicPanelOfNormals}.  See docs/CNVs/CNV-methods.pdf for details.
//     */
//    public static double hetLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
//                                          final AllelicCount count, final AlleleFractionIndicator indicator,
//                                          final AllelicPanelOfNormals allelicPoN) {
//        final SimpleInterval site = count.getInterval();
//        final double alpha = allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON) ? parameters.getAlpha() : allelicPoN.getAlpha(site);
//        final double beta = allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON) ? parameters.getBeta() : allelicPoN.getBeta(site);
//        final double pi = parameters.getOutlierProbability();
//        final int a = count.getAltReadCount();
//        final int r = count.getRefReadCount();
//        return hetLogLikelihood(alpha, beta, minorFraction, pi, indicator, a, r);
//    }
//
//    private static double hetLogLikelihood(final double alpha, final double beta, final double minorFraction,
//                                           final double outlierProbability, final AlleleFractionIndicator indicator,
//                                           final int a, final int r) {
//        if (indicator == AlleleFractionIndicator.OUTLIER) {
//            return log(outlierProbability) + LN_2;  //outlier distribution is flat in [0, 0.5]
//        } else {
//            final double fAlt = indicator == AlleleFractionIndicator.ALT_MINOR ? minorFraction : 1 - minorFraction;
//            return log((1 - outlierProbability) / 2) + logIntegralOverAllelicBias(alpha, beta, fAlt, a, r);
//        }
//    }
//
//    public static double collapsedHetLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
//                                                   final AllelicCount count, final AllelicPanelOfNormals allelicPoN) {
//        return GATKProtectedMathUtils.logSumExp(
//                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.ALT_MINOR, allelicPoN),
//                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.REF_MINOR, allelicPoN),
//                hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.OUTLIER, allelicPoN));
//    }
//
//    public static double segmentLogLikelihood(final AlleleFractionGlobalParameters parameters, final double minorFraction,
//                                              final Collection<AllelicCount> counts, final AllelicPanelOfNormals allelicPoN) {
//        return counts.stream().mapToDouble(c -> collapsedHetLogLikelihood(parameters, minorFraction, c, allelicPoN)).sum();
//    }
//
//    public static double logLikelihood(final AlleleFractionGlobalParameters parameters, final AlleleFractionState.MinorFractions minorFractions,
//                                       final AlleleFractionData data) {
//        return IntStream.range(0, data.getNumSegments())
//                .mapToDouble(s -> segmentLogLikelihood(parameters, minorFractions.get(s), data.getCountsInSegment(s), data.getPoN())).sum();
//    }
//
//    protected static double logIntegralOverAllelicBias(final double sampleDepth, final double sampleBias, final double fAlt, final int a, final int r) {
//
//    }
}