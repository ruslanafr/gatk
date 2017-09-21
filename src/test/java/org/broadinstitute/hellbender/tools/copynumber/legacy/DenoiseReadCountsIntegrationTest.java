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
    private static final File WGS_ANNOTATED_INTERVALS_FILE = new File(TEST_SUB_DIR, "denoise-read-counts-SM-74P4M-v1-chr20-downsampled-wgs-read-counts.readCounts.intervals.annotated.tsv");
    private static final File WGS_NO_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-no-gc.pon.hdf5");
    private static final File WGS_DO_GC_PON_FILE = new File(largeFileTestDir, "cnv_somatic_workflows_test_files/wgs-do-gc.pon.hdf5");

    @Test
    public void testWithNoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_TSV_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }

    @Test
    public void testHDF5InputWithNoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_HDF5_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }

    @Test
    public void testWithNoGCPonWithAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_TSV_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME, WGS_ANNOTATED_INTERVALS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }

    @Test
    public void testWithoutPoNWithAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_HDF5_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME, WGS_ANNOTATED_INTERVALS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names, intervals, and standardized and denoised copy ratios should be the same when no PoN is provided
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios, denoisedCopyRatios);
    }

    @Test
    public void testWithoutPoNWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_HDF5_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names, intervals, and standardized and denoised copy ratios should be the same when no PoN is provided
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios, denoisedCopyRatios);
    }

    @Test
    public void testWithDoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_TSV_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_DO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }

    @Test
    public void testHDF5InputWithDoGCPonWithoutAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_HDF5_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }

    @Test
    public void testWithDoGCPonWithAnnotations() {
        final File standardizedCRFile = createTempFile("test", ".standardizedCR.tsv");
        final File denoisedCRFile = createTempFile("test", ".denoisedCR.tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, WGS_READ_COUNTS_TSV_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_DO_GC_PON_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME, WGS_ANNOTATED_INTERVALS_FILE.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, standardizedCRFile.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, denoisedCRFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);

        //check that warning is thrown when both PoN with GC annotations and GC-annotated intervals are provided


        Assert.assertTrue(standardizedCRFile.exists());
        Assert.assertTrue(denoisedCRFile.exists());

        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(standardizedCRFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(denoisedCRFile);
        Assert.assertFalse(standardizedCopyRatios == denoisedCopyRatios);
        //sample names and intervals should be the same, standardized and denoised copy ratios should be different
        //we do not check correctness of the output standardized or denoised copy ratios
        Assert.assertEquals(standardizedCopyRatios.getSampleName(), denoisedCopyRatios.getSampleName());
        Assert.assertEquals(standardizedCopyRatios.getIntervals(), denoisedCopyRatios.getIntervals());
        Assert.assertNotEquals(standardizedCopyRatios.getCopyRatioValues(), denoisedCopyRatios.getCopyRatioValues());
    }
}