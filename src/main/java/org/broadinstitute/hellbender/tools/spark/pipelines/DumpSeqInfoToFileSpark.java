package org.broadinstitute.hellbender.tools.spark.pipelines;

import java.util.StringJoiner;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;


@CommandLineProgramProperties(summary = "Remove SEQ and QUAL information from the BAM's reads and dump the information into a compressed text file given by --output",
        oneLineSummary = "Save SEQ and QUAL info from the given BAM into a file", programGroup = SparkProgramGroup.class)
@DocumentedFeature
@BetaFeature
public class DumpSeqInfoToFileSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            //https://github.com/broadinstitute/hellbender/issues/929
            throw new UserException("PrintReadsSpark: Only coordinate-sorted files are currently supported");
        }

        final JavaRDD<GATKRead> reads = getReads();

        final JavaRDD<String> readSeqInfo = reads.map(this::returnReadSeqInfo).distinct();
        readSeqInfo.coalesce(1,true).saveAsTextFile("readSeqInfo.gz", org.apache.hadoop.io.compress.GzipCodec.class);

    }

    private String returnReadSeqInfo(GATKRead read) {
        final SAMRecord readSAM = read.convertToSAMRecord(getHeaderForReads());
        final StringJoiner sj = new StringJoiner("\t");
        sj.add(read.getName())
                .add(Integer.toString(readSAM.getFlags()))
                .add(read.getBasesString())
                .add(readSAM.getBaseQualityString());
        return sj.toString();
    }
}
