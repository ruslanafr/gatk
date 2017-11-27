package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;
import java.util.Arrays;

@CommandLineProgramProperties(summary = "Remove SEQ and QUAL information from the BAM's reads, this information can be found in the reference alignment BAM",
        oneLineSummary = "Edits reads to remove SEQ and QUAL, and prints them using Spark", programGroup = SparkProgramGroup.class)
@DocumentedFeature
@BetaFeature
public class RemoveSeqInfoSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second alignment BAM that will be used as reference", shortName = "RB", fullName = "refBAM", optional = false)
    public String refBAM;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String output;


    private static Tuple2<Tuple2<String, Integer>, GATKRead> getNamedReadTuple (GATKRead read, SAMFileHeader header){
        //extract read name
        final String readName = read.getName();
        //turn the read into SAMRecord
        final SAMRecord readSAM = read.convertToSAMRecord(header);
        //extract the flag value from the SAMRecord
        final Integer readFlag = readSAM.getFlags();
        //figure out if the read is first or second in template
        final Integer mate_direction_flag_options = 0x9f0;
        final Integer finalReadFlag = readFlag & mate_direction_flag_options;

        final Tuple2<String, Integer> nameTuple = new Tuple2<>(readName, finalReadFlag);

        return new Tuple2<>(nameTuple, read);
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            //https://github.com/broadinstitute/hellbender/issues/929
            throw new UserException("PrintReadsSpark: Only coordinate-sorted files are currently supported");
        }

        //extract reads from the reference file
        final ReadsSparkSource refReadsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        final SAMFileHeader refReadsHeader = refReadsSource.getHeader(refBAM, null);
        final JavaRDD<GATKRead> ref_reads =  refReadsSource.getParallelReads(refBAM, null, null, bamPartitionSplitSize);
        //extract reads from the bam we need to change
        final JavaRDD<GATKRead> reads = getReads();

        JavaPairRDD<Tuple2<String, Integer>, Iterable<GATKRead>> namedRefReads = ref_reads.mapToPair(read -> getNamedReadTuple(read, refReadsHeader)).groupByKey();

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> namedReadsToFix = reads.mapToPair(read -> getNamedReadTuple(read, getHeaderForReads()));

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> readsNotInReference = namedReadsToFix.subtractByKey(namedRefReads);

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> readsThatExistInReference = namedReadsToFix.subtractByKey(readsNotInReference);

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> transformedReads = readsThatExistInReference.mapValues(this::removeSeqData);

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> finalListOfReadsData = transformedReads.union(readsNotInReference);

        JavaRDD<GATKRead> finalListOfReads = finalListOfReadsData.map((Tuple2<Tuple2<String, Integer>, GATKRead> pair) -> pair._2);

        writeReads(ctx, output, finalListOfReads);
    }

    private GATKRead removeSeqData(GATKRead read) {
        if (read.isSupplementaryAlignment())
            return read;
        final byte[] asterisk = new byte[read.getBases().length];
        Arrays.fill(asterisk, "N".getBytes()[0]);
        read.setBases(asterisk);
        //read.setBaseQualities(asterisk);
        return read;
    }


}
