package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by tsato on 8/1/17.
 */
public class CollectDataForReadOrientationFilterIntegrationTest extends CommandLineProgramTest {

    /**
     * Test the tool on real bam to make sure that it does not crash
     */
    @Test
    public void testOnRealBam() {
        final File refTable = createTempFile("ref", ".table");
        final File altTable = createTempFile("alt", ".table");

        final String[] args = {
                "-R", v37_chr17_1Mb_Reference,
                "-I", NA12878_chr17_1k_BAM,
                "-alt_table", altTable.getAbsolutePath(),
                "-ref_table", refTable.getAbsolutePath()
        };

        runCommandLine(args);


    }

    @Test
    public void testOnSyntheticBam() throws IOException{
        final File refTable = createTempFile("ref", ".table");
        final File altTable = createTempFile("alt", ".table");

        final File samFile = createSyntheticSam();

        final String[] args = {
                "-R", v37_chr17_1Mb_Reference,
                "-I", samFile.getAbsolutePath(),
                "-alt_table", altTable.getAbsolutePath(),
                "-ref_table", refTable.getAbsolutePath()
        };

        double yo = 3;

    }

    private File createSyntheticSam() throws IOException {
        final List<Allele> alleles = Arrays.asList(Allele.create((byte) 'C', true), Allele.create((byte) 'A', false));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);

        // create a sam header
        final int numChromosomes = 1;
        final int startingChromosome = 1;
        final int chromosomeSize = 1000;
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(
                numChromosomes, startingChromosome, chromosomeSize);

        // create a sample list
        final int chromosomeIndex = 1;
        final String sampleName = "meganshand";
        final SampleList sampleList = new IndexedSampleList(sampleName);

        // specify chracteristics of reads

        // Reference bases at this position = CACTGTAAAA
        final int alignmentStart = 100_000_000;
        final int readLength = 10;
        final byte baseq = 20;
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseq);
        final int mapq = 60;
        final byte[] refReadBases = "CACTGTAAAA".getBytes();
        final byte[] altReadBases = "TACTGTAAAA".getBytes(); // C -> T at position 0
        final int numAltReads = 25;
        final int numRefReads = 75;
        final int numReads = numAltReads + numRefReads;

        // create reads
        final Map<String, List<GATKRead>> readMap = new LinkedHashMap<>();
        final List<GATKRead> reads = new ArrayList<>(numReads);


        // create ALT reads first
        for (int i = 0; i < numAltReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, chromosomeIndex,
                    alignmentStart, altReadBases, quals, "10M");
            read.setMappingQuality(mapq);
            reads.add(read);
        }

        // create REF reads
        for (int i = numAltReads; i < numReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, chromosomeIndex,
                    alignmentStart, refReadBases, quals, "10M");
            read.setMappingQuality(mapq);
            reads.add(read);
        }

        final File samFile = File.createTempFile("synthetic","sam");
        final SAMFileGATKReadWriter writer = new SAMFileGATKReadWriter(new SAMTextWriter(samFile));
        reads.forEach(r -> writer.addRead(r));

        return samFile;
    }

}