package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.broadinstitute.hellbender.tools.exome.TargetArgumentCollection;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Integration tests for {@link DenoiseReadCounts}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class DenoiseReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "copynumber/legacy/coverage";
    private static final File NORMAL_READ_COUNT_FILE = new File("/home/slee/working/ipython/wes.tsv");
    private static final String WGS_NO_GC_PON_FILE = largeFileTestDir + "cnv_somatic_workflows_test_files/wgs-no-gc.pon.hdf5";

    @Test
    public void testWGSWithNoGCPonWithoutAnnotations() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.tsv",
                "-" + LegacyCopyNumberArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME, WGS_NO_GC_PON_FILE,
                "-" + LegacyCopyNumberArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME, "wes_case.with-no-gc-pon_without-annot.ptn.tsv",
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "wes_case.with-no-gc-pon_without-annot.tn.tsv",
                "-" + LegacyCopyNumberArgument.NUMBER_OF_EIGENSAMPLES_SHORT_NAME, "10",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWESWithGCPonWithoutAnnotations() {
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
    public void testWESWithoutPoNWithAnnotations() {
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
    public void testWESWithoutPoNWithoutAnnotations() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.tsv",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_without-annot.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wes-pon-test/wes_case.without-pon_without-annot.tn.tsv",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGS5M() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test-5M/wgs_4.tsv",
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test-5M/wgs-5M.no-gc.pon",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test-5M/wgs_4.no-gc.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test-5M/wgs_4.no-gc.tn.tsv",
//                "-" + DenoiseReadCounts.NUMBER_OF_EIGENSAMPLES_SHORT_NAME, "10",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGS() {
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test/wgs_4.tsv",
                "-" + ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test/wgs.no-gc.pon",
                "-" + ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test/wgs_4.no-gc.ptn.tsv",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, "/home/slee/working/ipython/wgs-pon-test/wgs_4.no-gc.tn.tsv",
//                "-" + DenoiseReadCounts.NUMBER_OF_EIGENSAMPLES_SHORT_NAME, "10",
                "--" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }
}