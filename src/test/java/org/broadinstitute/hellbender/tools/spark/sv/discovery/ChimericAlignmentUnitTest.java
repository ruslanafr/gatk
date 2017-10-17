package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AlnModType;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple4;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider.*;

public class ChimericAlignmentUnitTest extends BaseTest {

    @Test(groups = "sv")
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = SVDiscoveryTestDataProvider.LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, AlnModType.NONE);

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test(groups = "sv")
    public void testFilterByNextAlignmentMayBeNovelInsertion() throws Exception {
        final AlignmentInterval overlappingRegion1 = new AlignmentInterval(new SimpleInterval("19", 48699881, 48700035), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0, 100, AlnModType.NONE);
        final AlignmentInterval overlappingRegion2 = new AlignmentInterval(new SimpleInterval("19", 48700584, 48700669), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0, 100, AlnModType.NONE);

        Assert.assertTrue(ChimericAlignment.nextAlignmentMayBeInsertion(overlappingRegion1, overlappingRegion2, 50, CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true));
    }

    @DataProvider(name = "forBooleanSeriesAndSerialization")
    private Object[][] createTestData() {
        final List<Object[]> data = new ArrayList<>(20);

        // simple inversion
        Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, true, false, false});

        testData = forSimpleInversionWithHom_leftPlus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.FORWARD_TO_REVERSE, true, false, false});

        testData = forSimpleInversionWithHom_leftMinus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.FORWARD_TO_REVERSE, false, false, false});

        testData = forSimpleInversionWithHom_rightPlus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, true, false, false});

        testData = forSimpleInversionWithHom_rightMinus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, false, false, false});

        // simple deletion
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // tandem duplication simple contraction
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // tandem duplication simple expansion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});

        // tandem duplication simple expansion with novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});



        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});


        // second test: contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});


        // third test: contraction from 3 units to 2 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});


        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false});


        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forBooleanSeriesAndSerialization", groups = "sv")
    public void test(final AlignmentInterval region1, final AlignmentInterval region2,
                     final StrandSwitch expectedStrandSwitch,
                     final boolean expectedIsForwardStrandRepresentation,
                     final boolean expectedIsLikelySimpleTranslocation,
                     final boolean expectedIsLikelyInvDup) {

        Assert.assertEquals(ChimericAlignment.determineStrandSwitch(region1, region2), expectedStrandSwitch);
        Assert.assertEquals(ChimericAlignment.isForwardStrandRepresentation(region1, region2, expectedStrandSwitch), expectedIsForwardStrandRepresentation);
        Assert.assertEquals(ChimericAlignment.isLikelySimpleTranslocation(region1, region2, expectedStrandSwitch), expectedIsLikelySimpleTranslocation);
        Assert.assertEquals(ChimericAlignment.isLikelyInvertedDuplication(region1, region2), expectedIsLikelyInvDup);

        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "dummyName");
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, chimericAlignment);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final ChimericAlignment roundTrip = (ChimericAlignment) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, chimericAlignment);
        Assert.assertEquals(roundTrip.hashCode(), chimericAlignment.hashCode());
    }
}