package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

public class AFCRSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPointsAlleleFraction;
    private final int numPointsCopyRatio;
    private final double meanMinorAlleleFraction;
    private final double meanLog2CopyRatio;

    public AFCRSegment(final SimpleInterval interval,
                       final int numPointsAlleleFraction,
                       final int numPointsCopyRatio,
                       final double meanMinorAlleleFraction,
                       final double meanLog2CopyRatio) {
        Utils.nonNull(interval);
        ParamUtils.isPositive(numPointsAlleleFraction, "Number of allele-fraction points must be positive.");
        ParamUtils.isPositive(numPointsCopyRatio, "Number of copy-ratio points must be positive.");
        this.interval = interval;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.meanMinorAlleleFraction = meanMinorAlleleFraction;
        this.meanLog2CopyRatio = meanLog2CopyRatio;
    }

    public AFCRSegment(final SimpleInterval interval,
                       final List<Double> alternateAlleleFractions,
                       final List<Double> denoisedLog2CopyRatios) {
        Utils.nonNull(interval);
        Utils.nonNull(alternateAlleleFractions);
        Utils.nonNull(denoisedLog2CopyRatios);
        this.interval = interval;
        numPointsAlleleFraction = alternateAlleleFractions.size();
        numPointsCopyRatio = denoisedLog2CopyRatios.size();
        meanMinorAlleleFraction = alternateAlleleFractions.stream().mapToDouble(aaf -> Math.min(aaf, 1. - aaf)).average().getAsDouble();
        meanLog2CopyRatio = denoisedLog2CopyRatios.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
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

    public int getNumPointsAlleleFraction() {
        return numPointsAlleleFraction;
    }

    public int getNumPointsCopyRatio() {
        return numPointsCopyRatio;
    }

    public double getMeanMinorAlleleFraction() {
        return meanMinorAlleleFraction;
    }

    public double getMeanLog2CopyRatio() {
        return meanLog2CopyRatio;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AFCRSegment that = (AFCRSegment) o;

        return numPointsAlleleFraction == that.numPointsAlleleFraction &&
                numPointsCopyRatio == that.numPointsCopyRatio &&
                Double.compare(that.meanMinorAlleleFraction, meanMinorAlleleFraction) == 0 &&
                Double.compare(that.meanLog2CopyRatio, meanLog2CopyRatio) == 0 && interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPointsAlleleFraction;
        result = 31 * result + numPointsCopyRatio;
        temp = Double.doubleToLongBits(meanMinorAlleleFraction);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(meanLog2CopyRatio);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
