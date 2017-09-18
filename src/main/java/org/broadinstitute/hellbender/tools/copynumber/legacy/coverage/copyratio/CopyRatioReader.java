package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

final class CopyRatioReader extends TableReader<CopyRatio> {
    private final File file;

    CopyRatioReader(final File file) throws IOException {
        super(file);
        this.file = file;
    }

    String getSampleName() {
        final List<String> sampleNameCommentLines = new ArrayList<>();
        try (final XReadLines reader = new XReadLines(file)) {
            for (final String line : reader) {
                if (!line.startsWith(CopyRatioCollection.SAMPLE_NAME_COMMENT_PREFIX)) {
                    break;
                }
                sampleNameCommentLines.add(line);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file);
        }
        if (sampleNameCommentLines.size() != 1) {
            throw new UserException.BadInput(String.format("File does not contain one sample name specified by %s.",
                    CopyRatioCollection.SAMPLE_NAME_COMMENT_PREFIX));
        }
        return sampleNameCommentLines.get(0).replace(CopyRatioCollection.SAMPLE_NAME_COMMENT_PREFIX, "");
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
