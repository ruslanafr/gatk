package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class AlleleFractionState extends ParameterizedState<AlleleFractionParameter> {
    public static final class MinorFractions {
        private final List<Double> minorFractions;

        public MinorFractions(final int numSegments) {
            minorFractions = new ArrayList<>(numSegments);
        }

        public MinorFractions(final List<Double> minorFractions) {
            Utils.nonNull(minorFractions);
            Utils.validateArg(minorFractions.stream().allMatch(f -> 0. <= f && f <= 0.5),
                    "List of minor fractions contains an element not in [0, 0.5].");
            this.minorFractions = new ArrayList<>(minorFractions);
        }

        public void add(final double minorFraction) {
            ParamUtils.inRange(minorFraction, 0., 0.5, "Minor fraction must be in [0, 0.5].");
            minorFractions.add(minorFraction);
        }

        public double get(final int index) {
            return minorFractions.get(index);
        }
    }

    public AlleleFractionState(final double sampleDepth,
                               final double sampleBias,
                               final double outlierProbability,
                               final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(AlleleFractionParameter.SAMPLE_DEPTH,
                        ParamUtils.isPositiveOrZero(sampleDepth, "Sample depth must be non-negative.")),
                new Parameter<>(AlleleFractionParameter.SAMPLE_BIAS,
                        ParamUtils.isPositiveOrZero(sampleBias, "Sample bias must be non-negative.")),
                new Parameter<>(AlleleFractionParameter.OUTLIER_PROBABILITY,
                        ParamUtils.inRange(outlierProbability, 0., 1., "Outlier probability must be in [0, 1].")),
                new Parameter<>(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractions)));
    }


    public double sampleDepth() {
        return get(AlleleFractionParameter.SAMPLE_DEPTH, Double.class);
    }

    public double sampleBias() {
        return get(AlleleFractionParameter.SAMPLE_BIAS, Double.class);
    }

    public double outlierProbability() {
        return get(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class);
    }

    public double segmentMinorFraction(final int segment) {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class).get(segment);
    }

    public AlleleFractionGlobalParameters globalParameters() {
        return new AlleleFractionGlobalParameters(sampleDepth(), sampleBias(), outlierProbability());
    }

    public MinorFractions minorFractions() {
        return get(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, MinorFractions.class);
    }
}