package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Integration test for {@link CallCopyRatioSegments}.
 */
public final class CallCopyRatioSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File(toolsTestDir, "copynumber/legacy/coverage/caller");
    private static final File TEST_DENOISED_COPY_RATIOS = new File(TEST_DIR, "call-copy-ratio-segments-denoised-copy-ratios.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR, "call-copy-ratio-segments-segments.tsv");

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, TEST_DENOISED_COPY_RATIOS.getAbsolutePath(),
                "-" + LegacyCopyNumberArgument.SEGMENTS_FILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments = new CalledCopyRatioSegmentCollection(outputFile);
        Assert.assertEquals(calledCopyRatioSegments.getRecords().stream().map(s -> s.getCall().getOutputString()).toArray(), new String[] {"+", "-", "0", "0"});
    }
}
