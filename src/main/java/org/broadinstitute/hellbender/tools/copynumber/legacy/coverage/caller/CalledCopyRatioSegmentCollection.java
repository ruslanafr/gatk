package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
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

public class CalledCopyRatioSegmentCollection {
    enum CalledCopyRatioSegmentTableColumn {
        SAMPLE,
        CONTIG,
        START,
        END,
        NUM_POINTS,
        MEAN_LOG2_COPY_RATIO,
        CALL;

        static final TableColumnCollection COLUMNS = new TableColumnCollection(
                SAMPLE, CONTIG, START, END, NUM_POINTS, MEAN_LOG2_COPY_RATIO, CALL);
    }

    private final List<CalledCopyRatioSegment> segments;

    public CalledCopyRatioSegmentCollection(final List<CalledCopyRatioSegment> segments) {
        this.segments = Utils.nonNull(segments);
    }

    public List<SimpleInterval> getIntervals() {
        return segments.stream().map(CalledCopyRatioSegment::getInterval).collect(Collectors.toList());
    }

    public List<CalledCopyRatioSegment> getSegments() {
        return Collections.unmodifiableList(segments);
    }

    public void write(final File file,
                      final String sampleName) {
        Utils.nonNull(file);
        Utils.nonNull(sampleName);
        try (final TableWriter<CalledCopyRatioSegment> writer =
                     TableUtils.writer(file, CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn.COLUMNS,
                             (segment, dataLine) ->
                                     dataLine.append(sampleName)
                                             .append(segment.getContig())
                                             .append(segment.getStart())
                                             .append(segment.getEnd())
                                             .append(segment.getNumPoints())
                                             .append(segment.getMeanLog2CopyRatio())
                                             .append(segment.getCall().getOutputString()))) {
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
        final CalledCopyRatioSegmentCollection that = (CalledCopyRatioSegmentCollection) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}