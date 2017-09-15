package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionKernelSegmenter {
    private static final Logger logger = LogManager.getLogger(AlleleFractionKernelSegmenter.class);

    //Gaussian kernel for a specified variance; if variance is zero, use a linear kernel
    private static final Function<Double, BiFunction<Double, Double, Double>> kernel =
            variance -> variance == 0.
                    ? (x, y) -> x * y
                    : (x, y) -> FastMath.exp(-(x - y) * (x - y) / variance);

    private final int numPointsTotal;
    private final Map<String, List<SimpleInterval>> intervalsPerChromosome;
    private final Map<String, List<Double>> alternateAlleleFractionsPerChromosome;    //in log2 space

    public AlleleFractionKernelSegmenter(final AllelicCountCollection allelicCounts) {
        Utils.nonNull(allelicCounts);
        numPointsTotal = allelicCounts.getAllelicCounts().size();
        intervalsPerChromosome = allelicCounts.getAllelicCounts().stream()
                .map(AllelicCount::getInterval)
                .collect(Collectors.groupingBy(
                        SimpleInterval::getContig,
                        LinkedHashMap::new,
                        Collectors.mapping(Function.identity(), Collectors.toList())));
        final double[] alternateAlleleFractions = allelicCounts.getAllelicCounts().stream()
                .mapToDouble(ac -> ac.getRefReadCount() + ac.getAltReadCount() == 0
                        ? 0.
                        : (double) ac.getAltReadCount() / (ac.getRefReadCount() + ac.getAltReadCount()))
                .toArray();
        alternateAlleleFractionsPerChromosome = IntStream.range(0, allelicCounts.getAllelicCounts().size()).boxed()
                .map(i -> new ImmutablePair<>(
                        allelicCounts.getAllelicCounts().get(i).getContig(),
                        alternateAlleleFractions[i]))
                .collect(Collectors.groupingBy(
                        Pair::getKey,
                        LinkedHashMap::new,
                        Collectors.mapping(Pair::getValue, Collectors.toList())));
    }

    public AlleleFractionSegmentationResult findSegmentation(final int maxNumChangepointsPerChromosome,
                                                             final double kernelVariance,
                                                             final int kernelApproximationDimension,
                                                             final List<Integer> windowSizes,
                                                             final double numChangepointsPenaltyLinearFactor,
                                                             final double numChangepointsPenaltyLogLinearFactor) {
        ParamUtils.isPositiveOrZero(maxNumChangepointsPerChromosome, "Maximum number of changepoints must be non-negative.");
        ParamUtils.isPositiveOrZero(kernelVariance, "Variance of Gaussian kernel must be non-negative (if zero, a linear kernel will be used).");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(new HashSet<>(windowSizes).size() == windowSizes.size(), "Window sizes must all be unique.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLinearFactor,
                "Linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLogLinearFactor,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");

        logger.info(String.format("Finding changepoints in %d data points and %d chromosomes...",
                numPointsTotal, alternateAlleleFractionsPerChromosome.size()));

        //loop over chromosomes, find changepoints, and create allele-fraction segments
        final List<AlleleFractionSegmentationResult.AlleleFractionSegment> segments = new ArrayList<>();
        for (final String chromosome : alternateAlleleFractionsPerChromosome.keySet()) {
            logger.info(String.format("Finding changepoints in chromosome %s...", chromosome));
            final List<Double> alternateAlleleFractionsInChromosome = alternateAlleleFractionsPerChromosome.get(chromosome);

            final List<Integer> changepoints = new KernelSegmenter<>(alternateAlleleFractionsInChromosome)
                .findChangepoints(maxNumChangepointsPerChromosome, kernel.apply(kernelVariance), kernelApproximationDimension,
                        windowSizes, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, true);

            if (!changepoints.contains(alternateAlleleFractionsInChromosome.size())) {
                changepoints.add(alternateAlleleFractionsInChromosome.size() - 1);
            }
            int previousChangepoint = -1;
            for (final int changepoint : changepoints) {
                final int start = intervalsPerChromosome.get(chromosome).get(previousChangepoint + 1).getStart();
                final int end = intervalsPerChromosome.get(chromosome).get(changepoint).getEnd();
                final List<Double> alternateAlleleFractionsInSegment = alternateAlleleFractionsInChromosome.subList(
                        previousChangepoint + 1, changepoint + 1);
                segments.add(new AlleleFractionSegmentationResult.AlleleFractionSegment(
                        new SimpleInterval(chromosome, start, end),
                        alternateAlleleFractionsInSegment));
                previousChangepoint = changepoint;
            }
        }
        logger.info(String.format("Found %d segments in %d chromosomes.", segments.size(), alternateAlleleFractionsPerChromosome.keySet().size()));
        return new AlleleFractionSegmentationResult(segments);
    }
}
