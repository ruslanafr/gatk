package org.broadinstitute.hellbender.tools.copynumber.legacy;

import com.google.common.collect.ImmutableSet;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionSegmentationResult;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentationResult;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model segmented copy ratio from denoised read counts.",
        oneLineSummary = "Model segmented copy ratio from denoised read counts",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ModelSegments extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //filename tags for output
    protected static final String FILTERED_ALLELIC_COUNTS_FILE_SUFFIX = ".filtered.ac.tsv";
    protected static final String SEGMENTS_FILE_SUFFIX = ".seg";
    protected static final String COPY_RATIO_SEGMENTS_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    protected static final String ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX = ".af" + SEGMENTS_FILE_SUFFIX;
    protected static final String UNION_SEGMENTS_FILE_SUFFIX = ".union" + SEGMENTS_FILE_SUFFIX;
    protected static final String BEGIN_FIT_FILE_TAG = ".model-begin";
    protected static final String FINAL_FIT_FILE_TAG = ".model-final";
    protected static final String CR_PARAMETER_FILE_SUFFIX = ".cr.param";
    protected static final String AF_PARAMETER_FILE_SUFFIX = ".af.param";

    protected static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maxNumSegmentsPerChromosome";
    protected static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME = "maxNumSegsPerChr";

    protected static final String MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME = "minTotalAlleleCount";
    protected static final String MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME = "minAC";

    protected static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernelVarianceCopyRatio";
    protected static final String KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME = "kernVarCR";

    protected static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernelVarianceAlleleFraction";
    protected static final String KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME = "kernVarAF";

    protected static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernelApproximationDimension";
    protected static final String KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME = "kernApproxDim";

    protected static final String WINDOW_SIZES_LONG_NAME = "windowSizes";
    protected static final String WINDOW_SIZES_SHORT_NAME = "winSizes";

    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME = "numChangepointsPenaltyFactorCopyRatio";
    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME = "numChangeptsPenCR";

    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME = "numChangepointsPenaltyFactorAlleleFraction";
    protected static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME = "numChangeptsPenAF";

    protected static final String SMALL_COPY_RATIO_SEGMENT_THRESHOLD_LONG_NAME = "smallCopyRatioSegmentThreshold";
    protected static final String SMALL_COPY_RATIO_SEGMENT_THRESHOLD_SHORT_NAME = "smallCRSegTh";

    protected static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    protected static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    protected static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    protected static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    protected static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    protected static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    protected static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    protected static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "intervalThresholdCopyRatio";
    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "intervalThresholdAlleleFraction";
    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME = "maxNumIterationsSimSeg";
    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME = "maxIterSim";

    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME = "numIterationsSimSegPerFit";
    protected static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME = "numIterSimPerFit";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_FULL_NAME,
            shortName = LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = "Input file containing allelic counts (output of CollectAllelicCounts).",
            fullName = LegacyCopyNumberArgument.ALLELIC_COUNTS_FILE_FULL_NAME,
            shortName = LegacyCopyNumberArgument.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputAllelicCountsFile;

    @Argument(
            doc = "Prefix for output files. Will also be used as the sample name in downstream plots." +
                    "(Note: if this is a file path or contains slashes (/), " +
                    "the string after the final slash will be used as the sample name in downstream plots.)",
            fullName =  LegacyCopyNumberArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = LegacyCopyNumberArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            shortName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int maxNumSegmentsPerChromosome = 100;

    @Argument(
            doc = "Minimum total count required to include site in allele-fraction segmentation.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME,
            shortName = MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int minTotalAlleleCount = 10;

    @Argument(
            doc = "Variance of Gaussian kernel for copy-ratio segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_COPY_RATIO_LONG_NAME,
            shortName = KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double kernelVarianceCopyRatio = 0.;

    @Argument(
            doc = "Variance of Gaussian kernel for allele-fraction segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME,
            shortName = KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double kernelVarianceAlleleFraction = 1.;

    @Argument(
            doc = "Dimension of kernel approximation.  A subsample containing this number of data points " +
                    "will be taken from the copy-ratio profile and used to construct the approximation for each chromosome.  " +
                    "If the total number of data points in a chromosome is greater " +
                    "than this number, then all data points in the chromosome will be used.  " +
                    "Time complexity scales quadratically and space complexity scales linearly with this parameter.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            shortName = KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each data point to be a changepoint will be calculated " +
                    "assuming that it demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZES_LONG_NAME,
            shortName = WINDOW_SIZES_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private List<Integer> windowSizes = new ArrayList<>(Arrays.asList(8, 16, 32, 64, 128, 256));

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for copy-ratio segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double numChangepointsPenaltyFactorCopyRatio = 1.;

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for allele-fraction segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double numChangepointsPenaltyFactorAlleleFraction = 1.;

    @Argument(
            doc = "Threshold for small copy-ratio segment merging. " +
                    "If a segment contains strictly less than this number of copy-ratio points, " +
                    "it is considered small and will be merged with an adjacent segment.",
            fullName = SMALL_COPY_RATIO_SEGMENT_THRESHOLD_LONG_NAME,
            shortName = SMALL_COPY_RATIO_SEGMENT_THRESHOLD_SHORT_NAME,
            optional = true
    )
    private int smallCopyRatioSegmentThreshold = 3;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    private int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    private int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    private int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    private int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for copy-ratio similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    private double intervalThresholdCopyRatio = 4.;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    private double intervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for similar-segment merging.",
            fullName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME,
            optional = true
    )
    private int maxNumSimilarSegmentMergingIterations = 10;

    @Argument(
            doc = "Number of similar-segment--merging iterations per MCMC model refit. " +
                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
                    "Setting this to 0 will completely disable model refitting between iterations.)",
            fullName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME,
            shortName = NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME,
            optional = true
    )
    protected int numSimilarSegmentMergingIterationsPerFit = 0;

    //initialize data/segment variables, some of which may be optional
    private CopyRatioCollection denoisedCopyRatios = null;
    private CopyRatioSegmentationResult copyRatioSegments = CopyRatioSegmentationResult.NO_SEGMENTS;
    private AllelicCountCollection filteredAllelicCounts = null;
    private AllelicCountCollection allelicCounts = null;
    private AlleleFractionSegmentationResult alleleFractionSegments = AlleleFractionSegmentationResult.NO_SEGMENTS;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        validateArguments();

        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        //the string after the final slash in the output prefix (which may be an absolute file path) will be used as the sample name
        final String sampleName = outputPrefix.substring(outputPrefix.lastIndexOf("/") + 1);

        if (inputDenoisedCopyRatiosFile != null) {
            readDenoisedCopyRatios();                         //TODO remove use of ReadCountCollection
            performCopyRatioSegmentation();
            writeCopyRatioSegments(sampleName);
        }
        
        if (inputAllelicCountsFile != null) {
            readAndFilterAllelicCounts();
            performAlleleFractionSegmentation();
            writeAlleleFractionSegments(sampleName);
        }

//        //TODO legacy code is used for modelling here---replace with new models and python-based inference
//        //make Genome from input copy ratio and allele counts; need to trivially convert new AllelicCount to legacy AllelicCount
//        final Genome genome = new Genome(denoisedCopyRatios,
//                allelicCounts.getAllelicCounts().stream().map(ac -> new AllelicCount(ac.getInterval(), ac.getRefReadCount(), ac.getAltReadCount())).collect(Collectors.toList()));
//
//        logger.info("Combining copy-ratio and allele-fraction segments...");
//        final List<SimpleInterval> unionedSegments = SegmentUtils.unionSegments(
//                copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals(), genome);
//        logger.info("Number of segments after segment union: " + unionedSegments.size());
//        final File unionedSegmentsFile = new File(outputPrefix + UNION_SEGMENTS_FILE_SUFFIX);
//        SegmentUtils.writeSegmentFileWithNumTargetsAndNumSNPs(unionedSegmentsFile, unionedSegments, genome);
//        logSegmentsFileWrittenMessage(unionedSegmentsFile);
//
//        logger.info(String.format("Merging small copy-ratio segments (threshold = %d)...", smallCopyRatioSegmentThreshold));
//        final SegmentedGenome segmentedGenomeWithSmallSegments = new SegmentedGenome(unionedSegments, genome);
//        final SegmentedGenome segmentedGenome = segmentedGenomeWithSmallSegments.mergeSmallSegments(smallCopyRatioSegmentThreshold);
//        logger.info("Number of segments after small-segment merging: " + segmentedGenome.getSegments().size());
//
//        //initial MCMC model fitting performed by ACNVModeller constructor
//        final ACNVModeller modeller = new ACNVModeller(segmentedGenome,
//                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);
//
//        //write initial segments and parameters to file
//        writeACNVModeledSegmentAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);
//
//        //similar-segment merging (segment files are output for each merge iteration)
//        performSimilarSegmentMergingStep(modeller);
//
//        //write final segments and parameters to file
//        writeACNVModeledSegmentAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);
//
//        ctx.setLogLevel(originalLogLevel);
//        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }

    private void validateArguments() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised copy-ratio profile file or an allelic-counts file.");
        if (inputDenoisedCopyRatiosFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
    }

    private void readDenoisedCopyRatios() {
        logger.info(String.format("Reading denoised copy-ratio profile file (%s)...", inputDenoisedCopyRatiosFile));
        denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
    }

    private void performCopyRatioSegmentation() {
        logger.info("Starting copy-ratio segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        copyRatioSegments = new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorCopyRatio, numChangepointsPenaltyFactorCopyRatio);
    }

    private void writeCopyRatioSegments(final String sampleName) {
        final File copyRatioSegmentsFile = new File(outputPrefix + COPY_RATIO_SEGMENTS_FILE_SUFFIX);
        copyRatioSegments.write(copyRatioSegmentsFile, sampleName);
        logSegmentsFileWrittenMessage(copyRatioSegmentsFile);
    }

    private void readAndFilterAllelicCounts() {
        logger.info(String.format("Reading allelic-counts file (%s)...", inputAllelicCountsFile));
        final AllelicCountCollection unfilteredAllelicCounts = new AllelicCountCollection(inputAllelicCountsFile);
        logger.info(String.format("Filtering allelic counts with total count less than %d...", minTotalAlleleCount));
        filteredAllelicCounts = new AllelicCountCollection(unfilteredAllelicCounts.getAllelicCounts().stream()
                .filter(ac -> ac.getRefReadCount() + ac.getAltReadCount() >= minTotalAlleleCount)
                .collect(Collectors.toList()));
        final File filteredAllelicCountsFile = new File(outputPrefix + FILTERED_ALLELIC_COUNTS_FILE_SUFFIX);
        filteredAllelicCounts.write(filteredAllelicCountsFile);
        logger.info(String.format("Retained %d / %d sites after filtering on total count...",
                filteredAllelicCounts.getAllelicCounts().size(), unfilteredAllelicCounts.getAllelicCounts().size()));
        logger.info(String.format("Filtered allelic counts written to %s.", filteredAllelicCountsFile));
    }

    private void performAlleleFractionSegmentation() {
        logger.info("Starting allele-fraction segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        alleleFractionSegments = new AlleleFractionKernelSegmenter(filteredAllelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorAlleleFraction, numChangepointsPenaltyFactorAlleleFraction);
    }

    private void writeAlleleFractionSegments(final String sampleName) {
        final File alleleFractionSegmentsFile = new File(outputPrefix + ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX);
        alleleFractionSegments.write(alleleFractionSegmentsFile, sampleName);
        logSegmentsFileWrittenMessage(alleleFractionSegmentsFile);
    }

//    private void performSimilarSegmentMergingStep(final ACNVModeller modeller) {
//        logger.info("Merging similar segments...");
//        logger.info("Initial number of segments before similar-segment merging: " + modeller.getACNVModeledSegments().size());
//        //perform iterations of similar-segment merging until all similar segments are merged
//        for (int numIterations = 1; numIterations <= maxNumSimilarSegmentMergingIterations; numIterations++) {
//            logger.info("Similar-segment merging iteration: " + numIterations);
//            final int prevNumSegments = modeller.getACNVModeledSegments().size();
//            if (numSimilarSegmentMergingIterationsPerFit > 0 && numIterations % numSimilarSegmentMergingIterationsPerFit == 0) {
//                //refit model after this merge iteration
//                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, true);
//            } else {
//                //do not refit model after this merge iteration (deciles will be unspecified)
//                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, false);
//            }
//            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
//                break;
//            }
//        }
//        if (!modeller.isModelFit()) {
//            //make sure final model is completely fit (i.e., deciles are specified)
//            modeller.fitModel();
//        }
//        logger.info("Final number of segments after similar-segment merging: " + modeller.getACNVModeledSegments().size());
//    }
//
//    //write modeled segments and global parameters to file
//    private void writeACNVModeledSegmentAndParameterFiles(final ACNVModeller modeller, final String fileTag) {
//        final File modeledSegmentsFile = new File(outputPrefix + fileTag + SEGMENTS_FILE_SUFFIX);
//        modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);
//        logSegmentsFileWrittenMessage(modeledSegmentsFile);
//        final File copyRatioParameterFile = new File(outputPrefix + fileTag + CR_PARAMETER_FILE_SUFFIX);
//        final File alleleFractionParameterFile = new File(outputPrefix + fileTag + AF_PARAMETER_FILE_SUFFIX);
//        modeller.writeACNVModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
//    }

    private void logSegmentsFileWrittenMessage(final File file) {
        logger.info("Segments written to file " + file);
    }
}