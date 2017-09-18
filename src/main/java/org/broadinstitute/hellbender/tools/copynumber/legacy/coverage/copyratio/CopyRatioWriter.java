package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

final class CopyRatioWriter extends TableWriter<CopyRatio> {

    CopyRatioWriter(final File file) throws IOException {
        super(file, CopyRatioTableColumn.COLUMNS);
    }

    @Override
    protected void composeLine(final CopyRatio record, final DataLine dataLine) {
        Utils.nonNull(record);
        Utils.nonNull(dataLine);
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getCopyRatioValue());
    }
}
