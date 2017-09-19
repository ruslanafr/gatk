package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import org.broadinstitute.hellbender.exceptions.UserException;
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
public final class AlleleFractionSegmentCollection {
    enum AlleleFractionSegmentTableColumn {
        SAMPLE,
        CONTIG,
        START,
        END,
        NUM_POINTS,
        MINOR_ALLELE_FRACTION;

        static final TableColumnCollection COLUMNS = new TableColumnCollection(
                SAMPLE, CONTIG, START, END, NUM_POINTS, MINOR_ALLELE_FRACTION);
    }

    public static final AlleleFractionSegmentCollection NO_SEGMENTS = new AlleleFractionSegmentCollection(Collections.emptyList());

    private final List<AlleleFractionSegment> segments;

    public AlleleFractionSegmentCollection(final List<AlleleFractionSegment> segments) {
        this.segments = Utils.nonNull(segments);
    }

    public List<SimpleInterval> getIntervals() {
        return segments.stream().map(AlleleFractionSegment::getInterval).collect(Collectors.toList());
    }

    public void write(final File file,
                      final String sampleName) {
        Utils.nonNull(file);
        Utils.nonNull(sampleName);
        try (final TableWriter<AlleleFractionSegment> writer =
                     TableUtils.writer(file, AlleleFractionSegmentTableColumn.COLUMNS,
                             (segment, dataLine) ->
                                     dataLine.append(sampleName)
                                             .append(segment.getContig())
                                             .append(segment.getStart())
                                             .append(segment.getEnd())
                                             .append(segment.getNumPoints())
                                             .append(segment.getMeanMinorAlleleFraction()))) {
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
        final AlleleFractionSegmentCollection that = (AlleleFractionSegmentCollection) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}