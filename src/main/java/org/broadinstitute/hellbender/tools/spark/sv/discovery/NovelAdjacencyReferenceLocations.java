package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.List;
import java.util.Objects;

/**
 * This class represents a pair of inferred genomic locations on the reference whose novel adjacency is generated
 * due to a simple SV event (in other words, a simple rearrangement between two genomic locations)
 * that is suggested by the input {@link AlignedContig},
 * and complications in pinning down the locations to exact base pair resolution.
 */
@DefaultSerializer(NovelAdjacencyReferenceLocations.Serializer.class)
public class NovelAdjacencyReferenceLocations {

    public final SimpleInterval leftJustifiedLeftRefLoc;
    public final SimpleInterval leftJustifiedRightRefLoc;

    public final StrandSwitch strandSwitch;
    public final BreakpointComplications complication;

    public NovelAdjacencyReferenceLocations(final ChimericAlignment chimericAlignment, final byte[] contigSequence,
                                            final SAMSequenceDictionary referenceDictionary) {

        // first get strand switch type, then get complications, finally use complications to justify breakpoints
        strandSwitch = chimericAlignment.strandSwitch;

        complication = new BreakpointComplications(chimericAlignment, contigSequence);

        final Tuple2<SimpleInterval, SimpleInterval> leftJustifiedBreakpoints = leftJustifyBreakpoints(chimericAlignment, complication, referenceDictionary);
        leftJustifiedLeftRefLoc = leftJustifiedBreakpoints._1();
        leftJustifiedRightRefLoc = leftJustifiedBreakpoints._2();
    }

    protected NovelAdjacencyReferenceLocations(final Kryo kryo, final Input input) {
        final String contig1 = input.readString();
        final int start1 = input.readInt();
        final int end1 = input.readInt();
        this.leftJustifiedLeftRefLoc = new SimpleInterval(contig1, start1, end1);
        final String contig2 = input.readString();
        final int start2 = input.readInt();
        final int end2 = input.readInt();
        this.leftJustifiedRightRefLoc = new SimpleInterval(contig2, start2, end2);

        this.strandSwitch = StrandSwitch.values()[input.readInt()];
        this.complication = kryo.readObject(input, BreakpointComplications.class);
    }

    private abstract static class BreakpointLeftAligner {

        protected String upstreamBreakpointRefContig;
        protected String dnstreamBreakpointRefContig;
        protected int upstreamBreakpointRefPos;
        protected int dnstreamBreakpointRefPos;

        final Tuple2<SimpleInterval, SimpleInterval> getLeftJustifiedBreakpoints() {
            return new Tuple2<>(new SimpleInterval(upstreamBreakpointRefContig, upstreamBreakpointRefPos, upstreamBreakpointRefPos),
                                new SimpleInterval(dnstreamBreakpointRefContig, dnstreamBreakpointRefPos, dnstreamBreakpointRefPos));
        }

        ///////////////
        abstract static class SameChrEventsBreakpointsAligner extends BreakpointLeftAligner {
            protected SameChrEventsBreakpointsAligner(final ChimericAlignment ca) {
                upstreamBreakpointRefContig = dnstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
            }
        }

        /**
         * For simple insertion/deletion, inversion, and dispersed duplications without explicit duplication annotation
         */
        final static class NoExplicitDupAnnotationBreakpointsAligner extends SameChrEventsBreakpointsAligner {
            NoExplicitDupAnnotationBreakpointsAligner(final ChimericAlignment ca, final BreakpointComplications complication) {

                super(ca);

                final int homologyLen = complication.getHomologyForwardStrandRep().length();
                final AlignmentInterval one = ca.regionWithLowerCoordOnContig,
                                        two = ca.regionWithHigherCoordOnContig;
                final SimpleInterval leftReferenceInterval, rightReferenceInterval;
                if (ca.isForwardStrandRepresentation) {
                    leftReferenceInterval  = one.referenceSpan;
                    rightReferenceInterval = two.referenceSpan;
                } else {
                    leftReferenceInterval  = two.referenceSpan;
                    rightReferenceInterval = one.referenceSpan;
                }
                if (ca.strandSwitch == StrandSwitch.NO_SWITCH) {

                    final List<AlignmentInterval> deOverlappedTempAlignments = ContigAlignmentsModifier.removeOverlap(one, two, null);
                    final AlignmentInterval deOverlappedOne = deOverlappedTempAlignments.get(0),
                                            deOverlappedTwo = deOverlappedTempAlignments.get(1);
                    final boolean isLikelyDispersedDuplication =
                            deOverlappedOne.referenceSpan.getStart() > deOverlappedTwo.referenceSpan.getStart() == deOverlappedOne.forwardStrand;
                    if (isLikelyDispersedDuplication) { // left and right seem to be flipped but that's the feature of dispersed duplication
                        upstreamBreakpointRefPos = rightReferenceInterval.getStart();
                        dnstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                    } else {
                        upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                        dnstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;
                    }

                } else if (ca.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE){
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                    dnstreamBreakpointRefPos = rightReferenceInterval.getEnd();
                } else {
                    upstreamBreakpointRefPos = leftReferenceInterval.getStart() - 1;
                    dnstreamBreakpointRefPos = rightReferenceInterval.getStart() + homologyLen - 1;
                }

            }
        }

        final static class TanDupBreakpointsAligner extends SameChrEventsBreakpointsAligner {

            TanDupBreakpointsAligner(final ChimericAlignment ca,
                                     final BreakpointComplications complication) {
                super(ca);

                final int homologyLen = complication.getHomologyForwardStrandRep().length();

                final SimpleInterval leftReferenceInterval, rightReferenceInterval;
                if (ca.isForwardStrandRepresentation) {
                    leftReferenceInterval  = ca.regionWithLowerCoordOnContig.referenceSpan;
                    rightReferenceInterval = ca.regionWithHigherCoordOnContig.referenceSpan;
                } else {
                    leftReferenceInterval  = ca.regionWithHigherCoordOnContig.referenceSpan;
                    rightReferenceInterval = ca.regionWithLowerCoordOnContig.referenceSpan;
                }
                if (complication.getDupSeqRepeatNumOnCtg() > complication.getDupSeqRepeatNumOnRef()) {
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen
                            - (complication.getDupSeqRepeatNumOnCtg() - complication.getDupSeqRepeatNumOnRef()) * complication.getDupSeqRepeatUnitRefSpan().size();
                } else {
                    upstreamBreakpointRefPos = leftReferenceInterval.getEnd() - homologyLen;
                }
                dnstreamBreakpointRefPos = rightReferenceInterval.getStart() - 1;;
            }
        }

        final static class InvDupBreakpointsAligner extends SameChrEventsBreakpointsAligner {
            InvDupBreakpointsAligner(final ChimericAlignment ca,
                                     final BreakpointComplications complication) {
                super(ca);
                upstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getStart() - 1;
                dnstreamBreakpointRefPos = complication.getDupSeqRepeatUnitRefSpan().getEnd();
            }
        }

        ///////////////

        /**
         * For computing exact and left adjusted breakpoint locations of inter-chromosome novel adjacency,
         * with or without strand switch.
         */
        final static class TransLocBreakpointsAligner extends BreakpointLeftAligner {

            TransLocBreakpointsAligner(final ChimericAlignment ca,
                                       final BreakpointComplications complication,
                                       final SAMSequenceDictionary referenceDictionary) {


                determineRefContigs(ca, referenceDictionary);

                extractRefPositions(ca, complication, referenceDictionary);
            }

            private void extractRefPositions(final ChimericAlignment ca, final BreakpointComplications complication,
                                             final SAMSequenceDictionary referenceDictionary) {
                final int homologyLen = complication.getHomologyForwardStrandRep().length();
                final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
                if (firstInPartner) {
                    switch (ca.strandSwitch) {
                        case NO_SWITCH:
                            if (ca.isForwardStrandRepresentation) {
                                upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                                dnstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                            } else {
                                upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                                dnstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            }
                            break;
                        case FORWARD_TO_REVERSE:
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            dnstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd();
                            break;
                        case REVERSE_TO_FORWARD:
                            upstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                            dnstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart() + homologyLen;
                            break;
                        default: throw new GATKException("Unseen strand switch case for: " + ca.onErrStringRep());
                    }
                } else {
                    switch (ca.strandSwitch) {
                        case NO_SWITCH:
                            if (ca.isForwardStrandRepresentation) {
                                upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                                dnstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            } else {
                                upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                                dnstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart();
                            }
                            break;
                        case FORWARD_TO_REVERSE:
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getEnd() - homologyLen;
                            dnstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getEnd();
                            break;
                        case REVERSE_TO_FORWARD:
                            upstreamBreakpointRefPos = ca.regionWithHigherCoordOnContig.referenceSpan.getStart();
                            dnstreamBreakpointRefPos = ca.regionWithLowerCoordOnContig.referenceSpan.getStart() + homologyLen;
                            break;
                        default: throw new GATKException("Unseen strand switch case for: " + ca.onErrStringRep());
                    }
                }
            }

            private void determineRefContigs(ChimericAlignment ca, SAMSequenceDictionary referenceDictionary) {
                final boolean firstInPartner = isFirstInPartner(ca, referenceDictionary);
                if (firstInPartner) {
                    upstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
                    dnstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
                } else {
                    upstreamBreakpointRefContig = ca.regionWithHigherCoordOnContig.referenceSpan.getContig();
                    dnstreamBreakpointRefContig = ca.regionWithLowerCoordOnContig.referenceSpan.getContig();
                }
            }

            private static boolean isFirstInPartner(final ChimericAlignment ca, final SAMSequenceDictionary referenceDictionary) {
                switch (ca.strandSwitch) {
                    case NO_SWITCH: return 0 > IntervalUtils.compareContigs(ca.regionWithLowerCoordOnContig.referenceSpan,
                            ca.regionWithHigherCoordOnContig.referenceSpan, referenceDictionary);
                    case FORWARD_TO_REVERSE: case REVERSE_TO_FORWARD:
                        return ca.isForwardStrandRepresentation;
                    default:
                        throw new GATKException("Unseen strand switch case for: " + ca.onErrStringRep());
                }
            }
        }
    }

    /**
     * Returns the reference coordinates of the left and right breakpoints implied by this chimeric alignment.
     * If there is homologous sequence represented in the alignments, it will be assigned to the side of the breakpoint
     * with higher reference coordinates.
     */
    @VisibleForTesting
    static Tuple2<SimpleInterval, SimpleInterval> leftJustifyBreakpoints(final ChimericAlignment ca,
                                                                         final BreakpointComplications complication,
                                                                         final SAMSequenceDictionary referenceDictionary) {

        final boolean sameChromosome = ca.regionWithLowerCoordOnContig.referenceSpan.getContig()
                .equals(ca.regionWithHigherCoordOnContig.referenceSpan.getContig());
        final Tuple2<SimpleInterval, SimpleInterval> leftAdjustedBreakpoints;
        if (sameChromosome) {
            if (complication.hasDuplicationAnnotation()) {
                if (complication.hasDupSeqButNoStrandSwitch()) { // tandem duplication
                    leftAdjustedBreakpoints = new BreakpointLeftAligner.TanDupBreakpointsAligner(ca, complication)
                            .getLeftJustifiedBreakpoints();
                } else { // inverted duplication
                    leftAdjustedBreakpoints = new BreakpointLeftAligner.InvDupBreakpointsAligner(ca, complication)
                            .getLeftJustifiedBreakpoints();
                }
            } else { // no explicit dup annotation
                leftAdjustedBreakpoints = new BreakpointLeftAligner.NoExplicitDupAnnotationBreakpointsAligner(ca, complication)
                        .getLeftJustifiedBreakpoints();
            }
        } else {
            leftAdjustedBreakpoints = new BreakpointLeftAligner.TransLocBreakpointsAligner(ca, complication, referenceDictionary)
                    .getLeftJustifiedBreakpoints();
        }

        validateInferredLocations(leftAdjustedBreakpoints._1, leftAdjustedBreakpoints._2, referenceDictionary, ca, complication);

        return leftAdjustedBreakpoints;
    }

    private static void validateInferredLocations(final SimpleInterval leftBreakpoint, final SimpleInterval rightBreakpoint,
                                                  final SAMSequenceDictionary referenceSequenceDictionary,
                                                  final ChimericAlignment ca, final BreakpointComplications complication) {

        Utils.validate(IntervalUtils.isBefore(leftBreakpoint, rightBreakpoint, referenceSequenceDictionary) ||
                        leftBreakpoint.equals(rightBreakpoint),
                "Inferred novel adjacency reference locations have left location after right location : left " +
                        leftBreakpoint + "\tright " + rightBreakpoint + "\t"  + ca.onErrStringRep() + "\n" + complication.toString());

        Utils.validate(leftBreakpoint.getEnd() <= referenceSequenceDictionary.getSequence(leftBreakpoint.getContig()).getSequenceLength(),
                "Inferred breakpoint beyond reference sequence length: inferred location: " + leftBreakpoint +
                        "\tref contig length: " + referenceSequenceDictionary.getSequence(leftBreakpoint.getContig()).getSequenceLength() + "\n"
                        + ca.onErrStringRep() + "\n" + complication.toString());
        Utils.validate(rightBreakpoint.getEnd() <= referenceSequenceDictionary.getSequence(rightBreakpoint.getContig()).getSequenceLength(),
                "Inferred breakpoint beyond reference sequence length: inferred location: " + rightBreakpoint +
                        "\tref contig length: " + referenceSequenceDictionary.getSequence(rightBreakpoint.getContig()).getSequenceLength() + "\n"
                        + ca.onErrStringRep() + "\n" + complication.toString());
    }

    @VisibleForTesting
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        final NovelAdjacencyReferenceLocations that = (NovelAdjacencyReferenceLocations) o;

        if (leftJustifiedLeftRefLoc != null ? !leftJustifiedLeftRefLoc.equals(that.leftJustifiedLeftRefLoc)
                : that.leftJustifiedLeftRefLoc != null)
            return false;
        if (leftJustifiedRightRefLoc != null ? !leftJustifiedRightRefLoc.equals(that.leftJustifiedRightRefLoc)
                : that.leftJustifiedRightRefLoc != null)
            return false;

        return strandSwitch.equals(that.strandSwitch) && complication.equals(that.complication);
    }

    @VisibleForTesting
    @Override
    public int hashCode() {
        return Objects.hash(leftJustifiedLeftRefLoc, leftJustifiedRightRefLoc, complication, 2659*strandSwitch.ordinal());
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeString(leftJustifiedLeftRefLoc.getContig());
        output.writeInt(leftJustifiedLeftRefLoc.getStart());
        output.writeInt(leftJustifiedLeftRefLoc.getEnd());
        output.writeString(leftJustifiedRightRefLoc.getContig());
        output.writeInt(leftJustifiedRightRefLoc.getStart());
        output.writeInt(leftJustifiedRightRefLoc.getEnd());
        output.writeInt(strandSwitch.ordinal());
        kryo.writeObject(output, complication);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<NovelAdjacencyReferenceLocations> {
        @Override
        public void write(final Kryo kryo, final Output output, final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public NovelAdjacencyReferenceLocations read(final Kryo kryo, final Input input, final Class<NovelAdjacencyReferenceLocations> klass ) {
            return new NovelAdjacencyReferenceLocations(kryo, input);
        }
    }

    /**
     * @return Intended for use in debugging and exception message only.
     */
    @Override
    public String toString() {
        return String.format("%s\t%s\t%s\t%s", leftJustifiedLeftRefLoc.toString(), leftJustifiedRightRefLoc.toString(),
                strandSwitch.name(), complication.toString());
    }
}
