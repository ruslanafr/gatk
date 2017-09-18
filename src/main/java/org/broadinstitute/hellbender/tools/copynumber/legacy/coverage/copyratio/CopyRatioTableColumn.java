package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

enum CopyRatioTableColumn {
    CONTIG,
    START,
    END,
    LOG2_COPY_RATIO;

    static final TableColumnCollection COLUMNS = new TableColumnCollection(
            CONTIG, START, END, LOG2_COPY_RATIO);
}
