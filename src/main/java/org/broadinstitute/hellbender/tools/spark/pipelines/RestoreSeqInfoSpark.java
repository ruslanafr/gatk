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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

/**
 * The goal of this program is to load two BAMs belonging to the same sample but aligned to two different references.
 * One of these BAMs does not contain any SEQ or QUAL info after running RemoveSeqInfoSpark on it.
 * This tool uses the second BAM to restore the missing SEQ and QUAL info.
 */
@CommandLineProgramProperties(summary = "Restores missing SEQ and QUAL info to the input BAM", oneLineSummary = "Restores missing SEQ and QUAL info to the input BAM",
        programGroup = SparkProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class RestoreSeqInfoSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the second alignment BAM that has the full information", shortName = "RB", fullName = "refBAM")
    public String refBAM;

    @Argument(doc="output BAM", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String output;


    private static Tuple2<Tuple2<String, Integer>, GATKRead> getNamedReadTuple (GATKRead read, SAMFileHeader header){
        final String readName = read.getName();
        final SAMRecord readSAM = read.convertToSAMRecord(header);
        final Integer readFlag = readSAM.getFlags();

        final Integer mate_direction_flag_options = 0x9f0;
        final Integer finalReadFlag = readFlag & mate_direction_flag_options;

        final Tuple2<String, Integer> nameTuple = new Tuple2<>(readName, finalReadFlag);

        return new Tuple2<>(nameTuple, read);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads_to_fix = getReads();

        final ReadsSparkSource refReadsSource = new ReadsSparkSource(ctx, readArguments.getReadValidationStringency());
        final SAMFileHeader refReadsHeader = refReadsSource.getHeader(refBAM, null);

        final JavaRDD<GATKRead> ref_reads =  refReadsSource.getParallelReads(refBAM, null, null, bamPartitionSplitSize);

        JavaPairRDD<Tuple2<String, Integer>, GATKRead> namedRefReads = ref_reads.mapToPair(read -> getNamedReadTuple(read, refReadsHeader))
                .reduceByKey((x,y) -> x);
        JavaPairRDD<Tuple2<String, Integer>, GATKRead> namedReadsToFix = reads_to_fix.mapToPair(read -> getNamedReadTuple(read, getHeaderForReads()));

        //put aside the reads that need to be in our final file but are not in the reference
        JavaPairRDD<Tuple2<String, Integer>, GATKRead> namedReadsNotInReference = namedReadsToFix.subtractByKey(namedRefReads);
        JavaRDD<GATKRead> readsNotInReference = namedReadsNotInReference.map((Tuple2<Tuple2<String, Integer>, GATKRead> pair) -> pair._2);


        JavaPairRDD<Tuple2<String, Integer>, GATKRead> namedReadsThatAreInReference = namedReadsToFix.subtractByKey(namedReadsNotInReference);

        JavaPairRDD<Tuple2<String, Integer>, Tuple2<GATKRead, GATKRead>> intersectionReads = namedReadsThatAreInReference.join(namedRefReads);

        JavaRDD<GATKRead> finalFixedReads = intersectionReads.mapValues((Tuple2<GATKRead, GATKRead> readPair) -> {
            //we found a match
            if (readPair._1.isSupplementaryAlignment()) {
                return readPair._1;
            }
            else {
                GATKRead refRead = readPair._2;
                GATKRead readToFix = readPair._1;

                readToFix.setBases(refRead.getBases());
                //readToFix.setBaseQualities(refRead.getBaseQualities());
                return readToFix;
            }
        }).map((Tuple2<Tuple2<String, Integer>, GATKRead> pair) -> pair._2);

        JavaRDD<GATKRead> finalReads = readsNotInReference.union(finalFixedReads);

        writeReads(ctx, output, finalReads);

    }

}


