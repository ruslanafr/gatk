package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the allele-fraction model containing the set of site alt and ref counts
 * and the grouping of sites into segments.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private final int numSegments;
    private final AllelicPanelOfNormals allelicPoN;
    private final List<AllelicCount> allelicCounts; // allelic counts indexed by site
    private final List<Integer> siteIndices; // the numbers 0, 1, 2 . . . N_sites
    private final List<Integer> startSitesPerSegment = new ArrayList<>();
    private final List<Integer> numSitesPerSegment = new ArrayList<>();

    public AlleleFractionData(final SegmentedGenome segmentedGenome) {
        this(segmentedGenome, AllelicPanelOfNormals.EMPTY_PON);
    }

    public AlleleFractionData(final SegmentedGenome segmentedGenome, final AllelicPanelOfNormals allelicPoN) {
        numSegments = segmentedGenome.getSegments().size();
        this.allelicPoN = allelicPoN;
        allelicCounts = new ArrayList<>();
        final List<SimpleInterval> segmentIntervals = segmentedGenome.getSegments();
//        final TargetCollection<AllelicCount> alleleCounts = segmentedGenome.getGenome().getSNPs();

        int startSite = 0;
        for (final SimpleInterval segment : segmentIntervals) {
            startSitesPerSegment.add(startSite);
//            final List<AllelicCount> countsInSegment = alleleCounts.targets(segment);
//            numSitesPerSegment.add(countsInSegment.size());
//            startSite += countsInSegment.size();
//            allelicCounts.addAll(countsInSegment);
        }

        siteIndices = IntStream.range(0, allelicCounts.size()).boxed().collect(Collectors.toList());
    }

    public AllelicPanelOfNormals getPoN() { return allelicPoN; }

    public List<AllelicCount> getAllelicCounts() { return Collections.unmodifiableList(allelicCounts); }

    public List<AllelicCount> getCountsInSegment(final int segment) {
        final int startInclusive = startSitesPerSegment.get(segment);
        final int endExclusive = startInclusive + numSitesPerSegment.get(segment);
        return Collections.unmodifiableList(allelicCounts.subList(startInclusive, endExclusive));
    }

    public List<Integer> getSitesInSegment(final int segment) {
        final int startInclusive = startSitesPerSegment.get(segment);
        final int endExclusive = startInclusive + numSitesPerSegment.get(segment);
        return Collections.unmodifiableList(siteIndices.subList(startInclusive, endExclusive));
    }

    public int getNumSitesInSegment(final int segment) {
        return numSitesPerSegment.get(segment);
    }

    public int getNumSegments() { return numSegments; }

    public AllelicCount getAllelicCount(final int site) { return allelicCounts.get(site); }

    public int getAltCount(final int site) { return allelicCounts.get(site).getAltReadCount(); }

    public int getRefCount(final int site) { return allelicCounts.get(site).getRefReadCount(); }

    public int getReadCount(final int site) { return getAltCount(site) + getRefCount(site); }
}
