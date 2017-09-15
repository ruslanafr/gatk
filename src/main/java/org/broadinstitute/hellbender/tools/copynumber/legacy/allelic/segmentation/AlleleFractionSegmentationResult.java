package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.SegmentTableColumn;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a legacy allele-fraction segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionSegmentationResult {
    //TODO update column headers; we keep the old ones for now
    private static final TableColumnCollection ALLELE_FRACTION_SEGMENT_FILE_TABLE_COLUMNS = new TableColumnCollection(
            SegmentTableColumn.SAMPLE,
            SegmentTableColumn.CONTIG,
            SegmentTableColumn.START,
            SegmentTableColumn.END,
            SegmentTableColumn.NUM_PROBES,
            SegmentTableColumn.MINOR_ALLELE_FRACTION);

    public static final AlleleFractionSegmentationResult NO_SEGMENTS = new AlleleFractionSegmentationResult(Collections.emptyList());

    private final List<AlleleFractionSegment> segments;

    static class AlleleFractionSegment {
        private final SimpleInterval interval;
        private final int numDataPoints;
        private final double meanMinorAlleleFraction;

        AlleleFractionSegment(final SimpleInterval interval,
                              final List<Double> alternateAlleleFractions) {
            Utils.nonNull(interval);
            Utils.nonEmpty(alternateAlleleFractions);
            this.interval = interval;
            numDataPoints = alternateAlleleFractions.size();
            meanMinorAlleleFraction = alternateAlleleFractions.stream().mapToDouble(aaf -> Math.min(aaf, 1. - aaf)).average().getAsDouble();
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
            if (numDataPoints != that.numDataPoints) {
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
            result = 31 * result + numDataPoints;
            temp = Double.doubleToLongBits(meanMinorAlleleFraction);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    public AlleleFractionSegmentationResult(final List<AlleleFractionSegment> segments) {
        this.segments = Collections.unmodifiableList(segments);
    }

    public List<SimpleInterval> getIntervals() {
        return segments.stream().map(s -> s.interval).collect(Collectors.toList());
    }

    public void write(final File file,
                      final String sampleName) {
        Utils.nonNull(file);
        Utils.nonNull(sampleName);
        try (final TableWriter<AlleleFractionSegment> writer =
                     TableUtils.writer(file, ALLELE_FRACTION_SEGMENT_FILE_TABLE_COLUMNS,
                             (segment, dataLine) ->
                                     dataLine.append(sampleName)
                                             .append(segment.interval.getContig())
                                             .append(segment.interval.getStart())
                                             .append(segment.interval.getEnd())
                                             .append(segment.numDataPoints)
                                             .append(segment.meanMinorAlleleFraction))) {
            writer.writeAllRecords(segments);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(file, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final AlleleFractionSegmentationResult that = (AlleleFractionSegmentationResult) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}