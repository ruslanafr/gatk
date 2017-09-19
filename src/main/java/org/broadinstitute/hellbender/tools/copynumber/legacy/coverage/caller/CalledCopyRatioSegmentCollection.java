package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

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

public class CalledCopyRatioSegmentCollection {

    private final String sampleName;
    private final List<CalledCopyRatioSegment> segments;

    public CalledCopyRatioSegmentCollection(final File inputFile) {
        IOUtils.canReadFile(inputFile);

        try (final CalledCopyRatioSegmentReader reader = new CalledCopyRatioSegmentReader(inputFile)) {
            sampleName = reader.getSampleName();
            segments = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public CalledCopyRatioSegmentCollection(final String sampleName,
                                            final List<CalledCopyRatioSegment> segments) {
        this.sampleName = Utils.nonNull(sampleName);
        this.segments = Utils.nonNull(segments);
    }

    public String getSampleName() {
        return sampleName;
    }

    public List<SimpleInterval> getIntervals() {
        return segments.stream().map(CalledCopyRatioSegment::getInterval).collect(Collectors.toList());
    }

    public List<CalledCopyRatioSegment> getSegments() {
        return Collections.unmodifiableList(segments);
    }

    public void write(final File outputFile) {
        try (final CalledCopyRatioSegmentWriter writer = new CalledCopyRatioSegmentWriter(outputFile, sampleName)) {
            writer.writeSampleName();
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
        final CalledCopyRatioSegmentCollection that = (CalledCopyRatioSegmentCollection) o;
        return segments.equals(that.segments);
    }

    @Override
    public int hashCode() {
        return segments.hashCode();
    }
}