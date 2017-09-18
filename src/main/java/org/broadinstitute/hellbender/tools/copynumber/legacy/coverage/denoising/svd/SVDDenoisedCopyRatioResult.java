package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Represents copy ratios for a sample that has been standardized and denoised by an {@link SVDReadCountPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SVDDenoisedCopyRatioResult {
    private final List<SimpleInterval> intervals;
    private final RealMatrix standardizedCopyRatioValues;
    private final RealMatrix denoisedCopyRatioValues;

    public SVDDenoisedCopyRatioResult(final List<SimpleInterval> intervals,
                                      final RealMatrix standardizedCopyRatioValues,
                                      final RealMatrix denoisedCopyRatioValues) {
        Utils.nonEmpty(intervals);
        Utils.nonNull(standardizedCopyRatioValues);
        Utils.nonNull(denoisedCopyRatioValues);
        Utils.validateArg(standardizedCopyRatioValues.getRowDimension() == 1,
                "Standardized copy-ratio values must contain only a single row.");
        Utils.validateArg(denoisedCopyRatioValues.getRowDimension() == 1,
                "Denoised copy-ratio values must contain only a single row.");
        Utils.validateArg(intervals.size() == standardizedCopyRatioValues.getColumnDimension(),
                "Number of intervals and columns in standardized copy-ratio values must match.");
        Utils.validateArg(intervals.size() == denoisedCopyRatioValues.getColumnDimension(),
                "Number of intervals and columns in denoised copy-ratio values must match.");
        this.intervals = intervals;
        this.standardizedCopyRatioValues = standardizedCopyRatioValues;
        this.denoisedCopyRatioValues = denoisedCopyRatioValues;
    }

    public void write(final File standardizedCopyRatiosFile,
                      final File denoisedCopyRatiosFile) {
        Utils.nonNull(standardizedCopyRatiosFile);
        Utils.nonNull(denoisedCopyRatiosFile);
        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(intervals, standardizedCopyRatioValues);
        standardizedCopyRatios.write(standardizedCopyRatiosFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(intervals, denoisedCopyRatioValues);
        denoisedCopyRatios.write(denoisedCopyRatiosFile);
    }
}
