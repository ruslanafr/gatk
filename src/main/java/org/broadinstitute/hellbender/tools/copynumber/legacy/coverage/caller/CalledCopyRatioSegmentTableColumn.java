package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

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
