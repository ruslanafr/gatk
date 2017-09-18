package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class CopyRatio implements Locatable {
    private final SimpleInterval interval;
    private final double copyRatioValue;

    public CopyRatio(final SimpleInterval interval,
                     final double copyRatioValue) {
        Utils.nonNull(interval);
        this.interval = interval;
        this.copyRatioValue = copyRatioValue;
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

    public double getCopyRatioValue() {
        return copyRatioValue;
    }
}
