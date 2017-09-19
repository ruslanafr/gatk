package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.NamedSampleFile;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

final class CopyRatioSegmentReader extends TableReader<CopyRatioSegment> implements NamedSampleFile {
    private final File file;

    CopyRatioSegmentReader(final File file) throws IOException {
        super(file);
        this.file = file;
    }

    String getSampleName() {
        return getSampleName(file);
    }

    @Override
    protected CopyRatioSegment createRecord(final DataLine dataLine) {
        Utils.nonNull(dataLine);
        try {
            final String contig = dataLine.get(CopyRatioSegmentTableColumn.CONTIG);
            final int start = dataLine.getInt(CopyRatioSegmentTableColumn.START);
            final int end = dataLine.getInt(CopyRatioSegmentTableColumn.END);
            final int numPoints = dataLine.getInt(CopyRatioSegmentTableColumn.NUM_POINTS);
            final double meanLog2CopyRatio = dataLine.getDouble(CopyRatioSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new CopyRatioSegment(interval, numPoints, meanLog2CopyRatio);
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("CopyRatioSegmentCollection file must have all columns specified.");
        }
    }
}
