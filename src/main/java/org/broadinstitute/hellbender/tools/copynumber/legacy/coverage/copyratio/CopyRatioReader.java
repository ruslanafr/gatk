package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

final class CopyRatioReader extends TableReader<CopyRatio> {

    CopyRatioReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected CopyRatio createRecord(final DataLine dataLine) {
        Utils.nonNull(dataLine);
        try {
            final String contig = dataLine.get(CopyRatioTableColumn.CONTIG);
            final int start = dataLine.getInt(CopyRatioTableColumn.START);
            final int end = dataLine.getInt(CopyRatioTableColumn.END);
            final double copyRatio = dataLine.getDouble(CopyRatioTableColumn.LOG2_COPY_RATIO);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new CopyRatio(interval, copyRatio);
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("CopyRatioCollection file must have all columns specified.");
        }
    }
}
