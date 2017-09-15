package org.broadinstitute.hellbender.tools.copynumber.legacy;

import com.google.common.collect.ImmutableSet;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentationResult;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
    //filename tags for output
    private static final String SEGMENTS_FILE_SUFFIX = ".seg";
    private static final String COPY_RATIO_SEGMENTS_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    private static final String ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX = ".af" + SEGMENTS_FILE_SUFFIX;
    private static final String UNION_SEGMENTS_FILE_SUFFIX = ".union" + SEGMENTS_FILE_SUFFIX;
    public static final String BEGIN_FIT_FILE_TAG = ".model-begin";
    public static final String FINAL_FIT_FILE_TAG = ".model-final";
    private static final String CR_PARAMETER_FILE_SUFFIX = ".cr.param";
    private static final String AF_PARAMETER_FILE_SUFFIX = ".af.param";

    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    private static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maxNumSegmentsPerChromosome";
    private static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME = "maxNumSegsPerChr";

    private static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernelVarianceCopyRatio";
    private static final String KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME = "kernVarCR";

    private static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernelVarianceAlleleFraction";
    private static final String KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME = "kernVarAF";

    private static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernelApproximationDimension";
    private static final String KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME = "kernApproxDim";

    private static final String WINDOW_SIZES_LONG_NAME = "windowSizes";
    private static final String WINDOW_SIZES_SHORT_NAME = "winSizes";

    private static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME = "numChangepointsPenaltyFactorCopyRatio";
    private static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME = "numChangeptsPenCR";

    private static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME = "numChangepointsPenaltyFactorAlleleFraction";
    private static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME = "numChangeptsPenAF";

    protected static final String SMALL_COPY_RATIO_SEGMENT_THRESHOLD_LONG_NAME = "smallCopyRatioSegmentThreshold";
    protected static final String SMALL_COPY_RATIO_SEGMENT_THRESHOLD_SHORT_NAME = "smallCRSegTh";

    private static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    private static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    private static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    private static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    private static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    private static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    private static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    private static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    private static final String INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "intervalThresholdCopyRatio";
    private static final String INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    private static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "intervalThresholdAlleleFraction";
    private static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    private static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME = "maxNumIterationsSimSeg";
    private static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME = "maxIterSim";

    private static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME = "numIterationsSimSegPerFit";
    private static final String NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_SHORT_NAME = "numIterSimPerFit";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = LegacyCopyNumberArgument.DENOISED_COPY_RATIO_PROFILE_FILE_FULL_NAME,
            shortName = LegacyCopyNumberArgument.DENOISED_COPY_RATIO_PROFILE_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatioProfileFile;

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
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
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
    protected int smallCopyRatioSegmentThreshold = 3;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for copy-ratio similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdCopyRatio = 4.;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for similar-segment merging.",
            fullName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumSimilarSegmentMergingIterations = 10;

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
    private ReadCountCollection denoisedCopyRatioProfile = null;
    private CopyRatioSegmentationResult copyRatioSegments = CopyRatioSegmentationResult.NO_SEGMENTS;
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

        if (inputDenoisedCopyRatioProfileFile != null) {
            readDenoisedCopyRatioProfile();
            performCopyRatioSegmentation();
            writeCopyRatioSegments(sampleName);
        } else {
            denoisedCopyRatioProfile = new ReadCountCollection(
                    Collections.emptyList(),
                    Collections.singletonList(sampleName),
                    new Array2DRowRealMatrix(0, 1));
        }
        
        if (inputAllelicCountsFile != null) {
            logger.info(String.format("Reading allelic-counts file (%s)...", inputAllelicCountsFile));
            allelicCounts = new AllelicCountCollection(inputAllelicCountsFile);
        }

        //TODO legacy code is used for modelling here---replace with new models and python-based inference
        //make Genome from input copy ratio and allele counts; need to trivially convert new AllelicCount to legacy AllelicCount
        final Genome genome = new Genome(denoisedCopyRatioProfile,
                allelicCounts.getCounts().stream().map(ac -> new AllelicCount(ac.getInterval(), ac.getRefReadCount(), ac.getAltReadCount())).collect(Collectors.toList()));

        logger.info("Combining copy-ratio and allele-fraction segments...");
        final List<SimpleInterval> unionedSegments = SegmentUtils.unionSegments(
                copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals(), genome);
        logger.info("Number of segments after segment union: " + unionedSegments.size());
        final File unionedSegmentsFile = new File(outputPrefix + UNION_SEGMENTS_FILE_SUFFIX);
        SegmentUtils.writeSegmentFileWithNumTargetsAndNumSNPs(unionedSegmentsFile, unionedSegments, genome);
        logSegmentsFileWrittenMessage(unionedSegmentsFile);

        logger.info("Merging small copy-ratio segments...");
        final SegmentedGenome segmentedGenomeWithSmallSegments = new SegmentedGenome(unionedSegments, genome);
        final SegmentedGenome segmentedGenome = segmentedGenomeWithSmallSegments.mergeSmallSegments(smallCopyRatioSegmentThreshold);
        logger.info("Number of segments after small-segment merging: " + segmentedGenome.getSegments().size());

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(segmentedGenome,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);

        //write initial segments and parameters to file
        writeACNVModeledSegmentAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);

        //similar-segment merging (segment files are output for each merge iteration)
        performSimilarSegmentMergingStep(modeller);

        //write final segments and parameters to file
        writeACNVModeledSegmentAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }

    private void validateArguments() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(inputDenoisedCopyRatioProfileFile == null && inputAllelicCountsFile == null,
                "Must provide at least a denoised copy-ratio profile file or an allelic-counts file.");
        if (inputDenoisedCopyRatioProfileFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatioProfileFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
    }

    private void readDenoisedCopyRatioProfile() {
        //TODO clean this up once updated ReadCountCollection is available, check sample names
        logger.info(String.format("Reading denoised copy-ratio profile file (%s)...", inputDenoisedCopyRatioProfileFile));
        try {
            denoisedCopyRatioProfile = ReadCountCollectionUtils.parse(inputDenoisedCopyRatioProfileFile);
        } catch (final IOException ex) {
            throw new UserException.BadInput("Could not read denoised copy-ratio profile file.");
        }
    }

    private void performCopyRatioSegmentation() {
        logger.info("Starting copy-ratio segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        copyRatioSegments = new CopyRatioKernelSegmenter(denoisedCopyRatioProfile)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorCopyRatio, numChangepointsPenaltyFactorCopyRatio);
    }

    private void writeCopyRatioSegments(final String sampleName) {
        final File copyRatioSegmentsFile = new File(outputPrefix + COPY_RATIO_SEGMENTS_FILE_SUFFIX);
        copyRatioSegments.write(copyRatioSegmentsFile, sampleName);
        logSegmentsFileWrittenMessage(copyRatioSegmentsFile);
    }

    private void performAlleleFractionSegmentation() {
        logger.info("Starting allele-fraction segmentation...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        alleleFractionSegments = new AlleleFractionKernelSegmenter(allelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorAlleleFraction, numChangepointsPenaltyFactorAlleleFraction);
    }

    private void writeAlleleFractionSegments(final String sampleName) {
        final File alleleFractionSegmentsFile = new File(outputPrefix + ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX);
        alleleFractionSegments.write(alleleFractionSegmentsFile, sampleName);
        logSegmentsFileWrittenMessage(alleleFractionSegmentsFile);
    }

    private void performSimilarSegmentMergingStep(final ACNVModeller modeller) {
        logger.info("Merging similar segments...");
        logger.info("Initial number of segments before similar-segment merging: " + modeller.getACNVModeledSegments().size());
        //perform iterations of similar-segment merging until all similar segments are merged
        for (int numIterations = 1; numIterations <= maxNumSimilarSegmentMergingIterations; numIterations++) {
            logger.info("Similar-segment merging iteration: " + numIterations);
            final int prevNumSegments = modeller.getACNVModeledSegments().size();
            if (numSimilarSegmentMergingIterationsPerFit > 0 && numIterations % numSimilarSegmentMergingIterationsPerFit == 0) {
                //refit model after this merge iteration
                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, true);
            } else {
                //do not refit model after this merge iteration (deciles will be unspecified)
                modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction, false);
            }
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
        }
        if (!modeller.isModelFit()) {
            //make sure final model is completely fit (i.e., deciles are specified)
            modeller.fitModel();
        }
        logger.info("Final number of segments after similar-segment merging: " + modeller.getACNVModeledSegments().size());
    }

    //write modeled segments and global parameters to file
    private void writeACNVModeledSegmentAndParameterFiles(final ACNVModeller modeller, final String fileTag) {
        final File modeledSegmentsFile = new File(outputPrefix + fileTag + SEGMENTS_FILE_SUFFIX);
        modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);
        logSegmentsFileWrittenMessage(modeledSegmentsFile);
        final File copyRatioParameterFile = new File(outputPrefix + fileTag + CR_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputPrefix + fileTag + AF_PARAMETER_FILE_SUFFIX);
        modeller.writeACNVModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void logSegmentsFileWrittenMessage(final File file) {
        logger.info("Segments written to file " + file);
    }
}