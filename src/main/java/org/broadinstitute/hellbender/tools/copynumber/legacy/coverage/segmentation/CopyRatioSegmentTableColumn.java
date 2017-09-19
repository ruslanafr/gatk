package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

enum CopyRatioSegmentTableColumn {
    SAMPLE,
    CONTIG,
    START,
    END,
    NUM_POINTS,
    MEAN_LOG2_COPY_RATIO;

    static final TableColumnCollection COLUMNS = new TableColumnCollection(
            SAMPLE, CONTIG, START, END, NUM_POINTS, MEAN_LOG2_COPY_RATIO);
}
