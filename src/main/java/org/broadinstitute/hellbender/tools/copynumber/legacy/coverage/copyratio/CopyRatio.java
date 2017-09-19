package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class CopyRatio implements Locatable {
    private final SimpleInterval interval;
    private final double log2CopyRatioValue;

    public CopyRatio(final SimpleInterval interval,
                     final double log2CopyRatioValue) {
        Utils.nonNull(interval);
        this.interval = interval;
        this.log2CopyRatioValue = log2CopyRatioValue;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public double getLog2CopyRatioValue() {
        return log2CopyRatioValue;
    }
}
