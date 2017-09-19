package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents a legacy copy-ratio segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioSegmentCollection {

    private final String sampleName;
    private final List<CopyRatioSegment> segments;


    public CopyRatioSegmentCollection(final File inputFile) {
        IOUtils.canReadFile(inputFile);

        try (final CopyRatioSegmentReader reader = new CopyRatioSegmentReader(inputFile)) {
            sampleName = reader.getSampleName();
            segments = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public CopyRatioSegmentCollection(final String sampleName,
                                      final List<CopyRatioSegment> segments) {
        this.sampleName = sampleName;
        this.segments = Utils.nonNull(segments);
    }

    public List<SimpleInterval> getIntervals() {
        return segments.stream().map(CopyRatioSegment::getInterval).collect(Collectors.toList());
    }

    public List<CopyRatioSegment> getSegments() {
        return Collections.unmodifiableList(segments);
    }

    public void write(final File outputFile) {
        try (final CopyRatioSegmentWriter writer = new CopyRatioSegmentWriter(outputFile, sampleName)) {
            writer.writeAllRecords(segments);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
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
        final CopyRatioSegmentCollection that = (CopyRatioSegmentCollection) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}