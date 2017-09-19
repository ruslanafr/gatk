package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public class AlleleFractionSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPoints;
    private final double meanMinorAlleleFraction;

    public AlleleFractionSegment(final SimpleInterval interval,
                                 final List<Double> alternateAlleleFractions) {
        Utils.nonNull(interval);
        Utils.nonEmpty(alternateAlleleFractions);
        this.interval = interval;
        numPoints = alternateAlleleFractions.size();
        meanMinorAlleleFraction = alternateAlleleFractions.stream().mapToDouble(aaf -> Math.min(aaf, 1. - aaf)).average().getAsDouble();
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

    public int getNumPoints() {
        return numPoints;
    }

    public double getMeanMinorAlleleFraction() {
        return meanMinorAlleleFraction;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final AlleleFractionSegment that = (AlleleFractionSegment) o;
        if (numPoints != that.numPoints) {
            return false;
        }
        if (Double.compare(that.meanMinorAlleleFraction, meanMinorAlleleFraction) != 0) {
            return false;
        }
        return interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPoints;
        temp = Double.doubleToLongBits(meanMinorAlleleFraction);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
