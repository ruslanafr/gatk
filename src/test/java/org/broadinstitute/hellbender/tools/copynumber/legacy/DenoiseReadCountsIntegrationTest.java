package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

/**
 * Integration tests for {@link DenoiseReadCounts}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class DenoiseReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/legacy/coverage/denoising";
    private static final File WGS_READ_COUNTS_TSV_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-SM-74P4M-v1-chr20-downsampled-wgs-read-counts.readCounts.tsv");
    private static final File WGS_READ_COUNTS_HDF5_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-SM-74P4M-v1-chr20-downsampled-wgs-read-counts.readCounts.hdf5");
    private static final File WGS_NO_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-no-gc.pon.hdf5");
    private static final File WGS_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-gc.pon.hdf5");

    @Test
    public void testWithNoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_TSV_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        //we do not check correctness of the output standardized or denoised copy ratios
    }

    @Test
    public void testHDF5InputWithNoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_HDF5_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        //we do not check correctness of the output standardized or denoised copy ratios
    }

    @Test
    public void testWithGCPonWithoutAnnotations() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.tsv",
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes.gc.pon",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.with-gc-pon_without-annot.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.with-gc-pon_without-annot.tn.tsv",
//                "-" + DenoiseReadCounts.NUMBER_OF_EIGENSAMPLES_SHORT_NAME, "10",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWithoutPoNWithAnnotations() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.tsv",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_with-annot.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_with-annot.tn.tsv",
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes.intervals.annot.tsv",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWithoutPoNWithoutAnnotations() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.tsv",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_without-annot.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_without-annot.tn.tsv",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }
}