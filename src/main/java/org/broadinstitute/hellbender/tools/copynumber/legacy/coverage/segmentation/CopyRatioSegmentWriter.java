package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

final class CopyRatioSegmentWriter extends TableWriter<CopyRatioSegment> {

    private final String sampleName;

    CopyRatioSegmentWriter(final File file,
                           final String sampleName) throws IOException {
        super(file, CopyRatioSegmentTableColumn.COLUMNS);
        this.sampleName = Utils.nonNull(sampleName);
    }

    @Override
    protected void composeLine(final CopyRatioSegment record, final DataLine dataLine) {
        Utils.nonNull(record);
        Utils.nonNull(dataLine);
        dataLine.append(sampleName)
                .append(record.getInterval().getContig())
                .append(record.getInterval().getStart())
                .append(record.getInterval().getEnd())
                .append(record.getNumPoints())
                .append(record.getMeanLog2CopyRatio());
    }
}
