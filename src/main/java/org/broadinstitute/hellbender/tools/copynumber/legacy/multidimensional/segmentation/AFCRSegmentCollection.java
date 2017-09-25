package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a legacy allele-fraction segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AFCRSegmentCollection extends TSVLocatableCollection<AFCRSegment> {
    enum AFCRSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_ALLELE_FRACTION,
        NUM_POINTS_COPY_RATIO,
        MEAN_MINOR_ALLELE_FRACTION,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, AFCRSegment> AFCR_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(AFCRSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(AFCRSegmentTableColumn.START);
        final int end = dataLine.getInt(AFCRSegmentTableColumn.END);
        final int numPointsAlleleFraction = dataLine.getInt(AFCRSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final int numPointsCopyRatio = dataLine.getInt(AFCRSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final double meanMinorAlleleFraction = dataLine.getDouble(AFCRSegmentTableColumn.MEAN_MINOR_ALLELE_FRACTION);
        final double meanLog2CopyRatio = dataLine.getDouble(AFCRSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new AFCRSegment(interval, numPointsAlleleFraction, numPointsCopyRatio, meanMinorAlleleFraction, meanLog2CopyRatio);
    };

    private static final BiConsumer<AFCRSegment, DataLine> AFCR_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getMeanMinorAlleleFraction())
                    .append(alleleFractionSegment.getMeanLog2CopyRatio());

    public AFCRSegmentCollection(final File inputFile) {
        super(inputFile, AFCRSegmentTableColumn.COLUMNS, AFCR_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, AFCR_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public AFCRSegmentCollection(final String sampleName,
                                 final List<AFCRSegment> AFCRSegments) {
        super(sampleName, AFCRSegments, AFCRSegmentTableColumn.COLUMNS, AFCR_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, AFCR_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}