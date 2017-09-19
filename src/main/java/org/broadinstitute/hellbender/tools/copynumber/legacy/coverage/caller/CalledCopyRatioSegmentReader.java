package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.NamedSampleFile;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

final class CalledCopyRatioSegmentReader extends TableReader<CalledCopyRatioSegment> implements NamedSampleFile {
    private final File file;

    CalledCopyRatioSegmentReader(final File file) throws IOException {
        super(file);
        this.file = file;
    }

    String getSampleName() {
        return getSampleName(file);
    }

    @Override
    protected CalledCopyRatioSegment createRecord(final DataLine dataLine) {
        Utils.nonNull(dataLine);
        try {
            final String contig = dataLine.get(CalledCopyRatioSegmentTableColumn.CONTIG);
            final int start = dataLine.getInt(CalledCopyRatioSegmentTableColumn.START);
            final int end = dataLine.getInt(CalledCopyRatioSegmentTableColumn.END);
            final int numPoints = dataLine.getInt(CalledCopyRatioSegmentTableColumn.NUM_POINTS);
            final double meanLog2CopyRatio = dataLine.getDouble(CalledCopyRatioSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
            final String callOutputString = dataLine.get(CalledCopyRatioSegmentTableColumn.CALL);
            final CalledCopyRatioSegment.Call call = Arrays.stream(CalledCopyRatioSegment.Call.values())
                    .filter(c -> c.getOutputString().equals(callOutputString)).findFirst().orElse(null);
            if (call == null) {
                throw new UserException.BadInput(String.format("Invalid call: %s", callOutputString));
            }
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new CalledCopyRatioSegment(new CopyRatioSegment(interval, numPoints, meanLog2CopyRatio), call);
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("CalledCopyRatioSegmentCollection file must have all columns specified.");
        }
    }
}
