package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd.HDF5SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.temporary.SimpleReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetAnnotation;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

/**
 * //TODO
 * Tool to create a panel of normals (PoN) given read counts for control samples.
 *
 * <p>
 * The input read counts are... 
 * </p>
 *
 * <h3>Examples</h3>
 *
 * <p>
 *     The following command is for either whole exome sequencing (WES) or whole genome sequencing (WGS) data.
 * </p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CreateReadCountPanelOfNormals \
 *   --input gc_corrected_coverages.tsv \
 *   --output panel_of_normals.pon
 * </pre>
 *
 * <p>
 * In addition to the resulting PoN, this command produces a .pon.removed_samples.txt file of samples removed for quality control (QC)
 * and a .pon.target_weights.txt file that gives the inverse variance per target that can optionally be passed to PerformSegmentation.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Create a panel of normals for read-count denoising given the read counts for samples in the panel.",
        oneLineSummary = "Create a panel of normals for read-count denoising.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class CreateReadCountPanelOfNormals extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //parameter names (some are public so they can be used in logging messages)
    public static final String MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME = "minimumIntervalMedianPercentile";
    static final String MINIMUM_INTERVAL_MEDIAN_PERCENTILE_SHORT_NAME = "minIntervalMedPct";
    public static final String MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME = "maximumZerosInSamplePercentage";
    static final String MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_SHORT_NAME = "maxZerosInSamplePct";
    public static final String MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME = "maximumZerosInIntervalPercentage";
    static final String MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_SHORT_NAME = "maxZerosInIntervalPct";
    public static final String EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME = "extremeSampleMedianPercentile";
    static final String EXTREME_SAMPLE_MEDIAN_PERCENTILE_SHORT_NAME = "extSampleMedPct";
    public static final String EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME = "extremeOutlierTruncationPercentile";
    static final String EXTREME_OUTLIER_TRUNCATION_PERCENTILE_SHORT_NAME = "extOutTruncPct";

    //default values for filtering (taken from ReCapSeg)
    private static final double DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE = 5.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE = 5.0;
    private static final double DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE = 5.0;
    private static final double DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE = 2.5;
    private static final double DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE = 0.1;

    private static final int DEFAULT_NUMBER_OF_EIGENSAMPLES = 20;

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples in the panel of normals.  " +
                    "Intervals must be identical and in the same order for all samples.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals (output of AnnotateTargets).  " +
                    "If provided, explicit GC correction will be performed before performing SVD.  " +
                    "Intervals must be identical to and in the same order as those in the input read-count files.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME,
            optional = true
    )
    private File annotatedIntervalsFile = null;

    @Argument(
            doc = "Output file name for the panel of normals.  " +
                    "Output is given in the HDF5 file format.  Contents can be viewed with the hdfview program.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputPanelOfNormalsFile;

    @Argument(
            doc = "Genomic intervals with a median (across samples) of fractional coverage (optionally corrected for GC bias) " +
                    "below this percentile are filtered out.  " +
                    "(This is the first filter applied.)",
            fullName = MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME,
            shortName = MINIMUM_INTERVAL_MEDIAN_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double minimumIntervalMedianPercentile = DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Samples with a fraction of zero-coverage genomic intervals above this percentage are filtered out.  " +
                    "(This is the second filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME,
            shortName = MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double maximumZerosInSamplePercentage = DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE;

    @Argument(
            doc = "Genomic intervals with a fraction of zero-coverage samples above this percentage are filtered out.  " +
                    "(This is the third filter applied.)",
            fullName = MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME,
            shortName = MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_SHORT_NAME,
            minValue = 0.,
            maxValue = 100.,
            optional = true
    )
    private double maximumZerosInIntervalPercentage = DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE;

    @Argument(
            doc = "Samples with a median (across genomic intervals) of fractional coverage normalized by genomic-interval medians  " +
                    "below this percentile or above the complementary percentile are filtered out.  " +
                    "(This is the fourth filter applied.)",
            fullName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME,
            shortName = EXTREME_SAMPLE_MEDIAN_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 50.,
            optional = true
    )
    private double extremeSampleMedianPercentile = DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE;

    @Argument(
            doc = "Fractional coverages normalized by genomic-interval medians that are " +
                    "below this percentile or above the complementary percentile are set to the corresponding percentile value.  " +
                    "(This is applied after all filters.)",
            fullName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME,
            shortName = EXTREME_OUTLIER_TRUNCATION_PERCENTILE_SHORT_NAME,
            minValue = 0.,
            maxValue = 50.,
            optional = true
    )
    private double extremeOutlierTruncationPercentile = DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE;

    @Argument(
            doc = "Number of eigensamples to use for truncated SVD and to store in the panel of normals.  " +
                    "The number of samples retained after filtering will be used instead if it is smaller than this.",
            fullName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            shortName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private int numEigensamplesRequested = DEFAULT_NUMBER_OF_EIGENSAMPLES;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        //validate parameters and parse optional parameters
        validateArguments();

        //get sample filenames
        final List<String> sampleFilenames = inputReadCountFiles.stream().map(File::getAbsolutePath).collect(Collectors.toList());

        //get intervals from the first read-count file to use as the canonical list of intervals
        //(this file is read again below, which is slightly inefficient but is probably not worth the extra code)
        final List<SimpleInterval> intervals = getIntervalsFromFirstReadCountFile(logger, inputReadCountFiles);

        //get GC content (null if not provided)
        final double[] intervalGCContent = validateIntervalGCContent(logger, intervals, annotatedIntervalsFile);

        //validate input read-count files (i.e., check intervals and that only integer counts are contained)
        //and aggregate as a RealMatrix with dimensions numIntervals x numSamples
        final RealMatrix readCountMatrix = constructReadCountMatrix(logger, inputReadCountFiles, intervals);

        //create the PoN
        logger.info("Creating the panel of normals...");
        HDF5SVDReadCountPanelOfNormals.create(outputPanelOfNormalsFile, getCommandLine(),
                readCountMatrix, sampleFilenames, intervals, intervalGCContent,
                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                extremeSampleMedianPercentile, extremeOutlierTruncationPercentile, numEigensamplesRequested, ctx);

        logger.info("Panel of normals successfully created.");
    }

    private void validateArguments() {
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        if (numEigensamplesRequested > inputReadCountFiles.size()) {
            logger.warn(String.format("Number of eigensamples (%d) is greater than the number of input samples (%d); " +
                            "the number of samples retained after filtering will be used instead.",
                    numEigensamplesRequested, inputReadCountFiles.size()));
        }
    }

    //TODO move GC-bias correction classes into copynumber package, clean up use of TargetCollection, and move this method into appropriate class
    //code is duplicated in DenoiseReadCounts for now
    private static double[] validateIntervalGCContent(final Logger logger,
                                                      final List<SimpleInterval> intervals,
                                                      final File annotatedIntervalsFile) {
        if (annotatedIntervalsFile == null) {
            logger.info("No GC-content annotations for intervals found; GC-bias correction will not be performed...");
            return null;
        }
        logger.info("Reading and validating GC-content annotations for intervals...");
        final TargetCollection<Target> annotatedIntervals = TargetArgumentCollection.readTargetCollection(annotatedIntervalsFile);
        Utils.validateArg(annotatedIntervals.targets().stream().map(Target::getInterval).collect(Collectors.toList()).equals(intervals),
                "Annotated intervals do not match intervals from first read-count file.");
        if (!annotatedIntervals.targets().stream().allMatch(t -> t.getAnnotations().hasAnnotation(TargetAnnotation.GC_CONTENT))) {
            throw new UserException.BadInput("At least one interval is missing a GC-content annotation.");
        }
        return annotatedIntervals.targets().stream().mapToDouble(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT)).toArray();
    }

    private static List<SimpleInterval> getIntervalsFromFirstReadCountFile(final Logger logger,
                                                                           final List<File> inputReadCountFiles) {
        final File firstReadCountFile = inputReadCountFiles.get(0);
        logger.info(String.format("Retrieving intervals from first read-count file (%s)...", firstReadCountFile));
        final SimpleReadCountCollection readCounts = parseReadCountFile(firstReadCountFile);
        return readCounts.getIntervals();
    }

    private RealMatrix constructReadCountMatrix(final Logger logger,
                                                final List<File> inputReadCountFiles,
                                                final List<SimpleInterval> intervals) {
        logger.info("Validating and aggregating input read-count files...");
        final int numSamples = inputReadCountFiles.size();
        final int numIntervals = intervals.size();
        final RealMatrix readCountMatrix = new Array2DRowRealMatrix(numSamples, numIntervals);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-count file %s (%d / %d)", inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleReadCountCollection readCounts = parseReadCountFile(inputReadCountFile);
            Utils.validateArg(readCounts.getIntervals().equals(intervals),
                    String.format("Intervals for read-count file %s do not match those in other read-count files.", inputReadCountFile));
            readCountMatrix.setRow(sampleIndex, readCounts.getReadCounts().getRow(0));
        }
        return readCountMatrix;
    }

    private static SimpleReadCountCollection parseReadCountFile(final File inputReadCountFile) {
        try (final HDF5File hdf5ReadCountFile = new HDF5File(inputReadCountFile)) {
            return SimpleReadCountCollection.read(hdf5ReadCountFile);
        } catch (final HDF5LibException e) {
            return SimpleReadCountCollection.read(inputReadCountFile);
        }
    }
}
