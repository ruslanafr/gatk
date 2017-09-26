package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

/**
 * TODO
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionModeller {
//    public static final double MAX_REASONABLE_SAMPLE_DEPTH = AlleleFractionInitializer.MAX_REASONABLE_SAMPLE_DEPTH;
//    public static final double MAX_REASONABLE_SAMPLE_BIAS = AlleleFractionInitializer.MAX_REASONABLE_SAMPLE_BIAS;
//
//    private final SegmentedGenome segmentedGenome;
//    private final ParameterizedModel<AlleleFractionParameter, AlleleFractionState, AlleleFractionData> model;
//    private final List<Double> sampleDepthSamples = new ArrayList<>();
//    private final List<Double> sampleBiasSamples = new ArrayList<>();
//    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
//    private final List<AlleleFractionState.MinorFractions> minorFractionsSamples = new ArrayList<>();
//    private final int numSegments;
//
//    public AlleleFractionModeller(final SegmentedGenome segmentedGenome, final AllelicPanelOfNormals allelicPoN) {
//        this.segmentedGenome = segmentedGenome;
//        final AlleleFractionData data = new AlleleFractionData(segmentedGenome, allelicPoN);
//        numSegments = data.getNumSegments();
//        final AlleleFractionState initialState = new AlleleFractionInitializer(data).getInitializedState();
//
//        // Initialization got us to the mode of the likelihood
//        // if we approximate conditionals as normal we can guess the width from the curvature at the mode and use as the slice-sampling widths
//        final AlleleFractionGlobalParameters initialParameters = initialState.globalParameters();
//        final AlleleFractionState.MinorFractions initialMinorFractions = initialState.minorFractions();
//
//        final double sampleDepthSamplingWidths = approximatePosteriorWidthAtMode(sampleDepth ->
//                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewSampleDepth(sampleDepth), initialMinorFractions, data), initialParameters.getSampleDepth());
//        final double sampleBiasSamplingWidths = approximatePosteriorWidthAtMode(sampleBias ->
//                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewSampleBias(sampleBias), initialMinorFractions, data), initialParameters.getSampleBias());
//        final double outlierProbabilitySamplingWidths = approximatePosteriorWidthAtMode(outlierProbability ->
//                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewOutlierProbability(outlierProbability), initialMinorFractions, data), initialParameters.getOutlierProbability());
//
//        final List<Double> minorFractionsSliceSamplingWidths = IntStream.range(0, numSegments).mapToDouble(segment ->
//                approximatePosteriorWidthAtMode(f -> AlleleFractionLikelihoods.segmentLogLikelihood(initialParameters, f, data.getCountsInSegment(segment), allelicPoN), initialMinorFractions.get(segment)))
//                .boxed().collect(Collectors.toList());
//
//        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> sampleDepthSampler =
//                new AlleleFractionSamplers.SampleDepthSampler(MAX_REASONABLE_SAMPLE_DEPTH, sampleDepthSamplingWidths);
//        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> sampleBiasSampler =
//                new AlleleFractionSamplers.SampleBiasSampler(MAX_REASONABLE_SAMPLE_BIAS, sampleBiasSamplingWidths);
//        final ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> outlierProbabilitySampler =
//                new AlleleFractionSamplers.OutlierProbabilitySampler(outlierProbabilitySamplingWidths);
//        final ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> minorFractionsSampler =
//                new AlleleFractionSamplers.MinorFractionsSampler(minorFractionsSliceSamplingWidths);
//
//        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
//                .addParameterSampler(AlleleFractionParameter.SAMPLE_DEPTH, sampleDepthSampler, Double.class)
//                .addParameterSampler(AlleleFractionParameter.SAMPLE_BIAS, sampleBiasSampler, Double.class)
//                .addParameterSampler(AlleleFractionParameter.OUTLIER_PROBABILITY, outlierProbabilitySampler, Double.class)
//                .addParameterSampler(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, minorFractionsSampler, AlleleFractionState.MinorFractions.class)
//                .build();
//    }
//
//    /**
//     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
//     * Gibbs sampling) to the collections held internally.  The current {@link AlleleFractionState} held internally is used
//     * to initialize the Markov Chain.
//     * @param numSamples    total number of samples per posterior
//     * @param numBurnIn     number of burn-in samples to discard
//     */
//    public void fitMCMC(final int numSamples, final int numBurnIn) {
//        //run MCMC
//        final GibbsSampler<AlleleFractionParameter, AlleleFractionState, AlleleFractionData> gibbsSampler = new GibbsSampler<>(numSamples, model);
//        gibbsSampler.runMCMC();
//
//        //update posterior samples
//        sampleDepthSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.SAMPLE_DEPTH, Double.class, numBurnIn));
//        sampleBiasSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.SAMPLE_BIAS, Double.class, numBurnIn));
//        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.OUTLIER_PROBABILITY, Double.class, numBurnIn));
//        minorFractionsSamples.addAll(gibbsSampler.getSamples(AlleleFractionParameter.MINOR_ALLELE_FRACTIONS, AlleleFractionState.MinorFractions.class, numBurnIn));
//    }
//
//    public List<Double> getsampleDepthSamples() {
//        return Collections.unmodifiableList(sampleDepthSamples);
//    }
//
//    public List<Double> getSampleBiasSamples() {
//        return Collections.unmodifiableList(sampleBiasSamples);
//    }
//
//    public List<Double> getOutlierProbabilitySamples() {
//        return Collections.unmodifiableList(outlierProbabilitySamples);
//    }
//
//    public List<AlleleFractionState.MinorFractions> getMinorFractionsSamples() {
//        return Collections.unmodifiableList(minorFractionsSamples);
//    }
//
//    public List<List<Double>> getMinorFractionSamplesBySegment() {
//        final List<List<Double>> result = new ArrayList<>();
//        for (int segment = 0; segment < numSegments; segment++) {
//            final List<Double> thisSegment = new ArrayList<>();
//            for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
//                thisSegment.add(sample.get(segment));
//            }
//            result.add(thisSegment);
//        }
//        return result;
//    }
//
//    /**
//     * Returns a list of {@link PosteriorSummary} elements summarizing the minor-allele-fraction posterior for each segment.
//     * Should only be called after {@link AlleleFractionModeller#fitMCMC(int, int)} has been called.
//     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
//     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
//     * @return                      list of {@link PosteriorSummary} elements summarizing the
//     *                              minor-allele-fraction posterior for each segment
//     */
//    public List<PosteriorSummary> getMinorAlleleFractionsPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
//        final int numSegments = segmentedGenome.getSegments().size();
//        final List<PosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
//        for (int segment = 0; segment < numSegments; segment++) {
//            final int j = segment;
//            final List<Double> minorFractionSamples =
//                    minorFractionsSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
//            posteriorSummaries.add(PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(minorFractionSamples, credibleIntervalAlpha, ctx));
//        }
//        return posteriorSummaries;
//    }
//
//    /**
//     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
//     * Should only be called after {@link AlleleFractionModeller#fitMCMC(int, int)} has been called.
//     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
//     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
//     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
//     */
//    public Map<AlleleFractionParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
//        final Map<AlleleFractionParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
//        posteriorSummaries.put(AlleleFractionParameter.SAMPLE_DEPTH, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(sampleDepthSamples, credibleIntervalAlpha, ctx));
//        posteriorSummaries.put(AlleleFractionParameter.SAMPLE_BIAS, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(sampleBiasSamples, credibleIntervalAlpha, ctx));
//        posteriorSummaries.put(AlleleFractionParameter.OUTLIER_PROBABILITY, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(outlierProbabilitySamples, credibleIntervalAlpha, ctx));
//        return posteriorSummaries;
//    }
//
//    //use width of a probability distribution given the position of its mode (estimated from Gaussian approximation) as step size
//    private static double approximatePosteriorWidthAtMode(final Function<Double, Double> logPDF, final double mode) {
//        final double absMode = Math.abs(mode);
//        final double epsilon = Math.min(1e-6, absMode / 2);    //adjust scale if mode is very near zero
//        final double defaultWidth = absMode / 10;              //if "mode" is not close to true mode of logPDF, approximation may not apply; just use 1/10 of absMode in this case
//        final double secondDerivative = (logPDF.apply(mode + epsilon) - 2 * logPDF.apply(mode) + logPDF.apply(mode - epsilon)) / (epsilon * epsilon);
//        return secondDerivative < 0 ? Math.sqrt(-1.0 / secondDerivative) : defaultWidth;
//    }
}