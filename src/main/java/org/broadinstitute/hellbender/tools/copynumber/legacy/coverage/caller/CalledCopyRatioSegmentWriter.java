package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.NamedSampleFile;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

public class CalledCopyRatioSegmentWriter extends TableWriter<CalledCopyRatioSegment> {
    private final String sampleName;

    CalledCopyRatioSegmentWriter(final File file,
                                 final String sampleName) throws IOException {
        super(file, CalledCopyRatioSegmentTableColumn.COLUMNS);
        this.sampleName = Utils.nonNull(sampleName);
    }

    void writeSampleName() {
        try {
            writeComment(NamedSampleFile.SAMPLE_NAME_COMMENT_PREFIX + sampleName);
        } catch (final IOException e) {
            throw new UserException("Could not write sample name.");
        }
    }

    @Override
    protected void composeLine(final CalledCopyRatioSegment record, final DataLine dataLine) {
        Utils.nonNull(record);
        Utils.nonNull(dataLine);
        dataLine.append(sampleName)
                .append(record.getInterval().getContig())
                .append(record.getInterval().getStart())
                .append(record.getInterval().getEnd())
                .append(record.getNumPoints())
                .append(record.getMeanLog2CopyRatio())
                .append(record.getCall().getOutputString());
    }
}