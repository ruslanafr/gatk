package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.EnumMap;
import java.util.List;
import java.util.stream.IntStream;

// statistics relevant to computing contamination and blacklisting possible LoH regions
class ContaminationStats {

    private double homAltCount;
    private double hetCount;
    private double expectedHomAltCount;
    private double expectedHetCount;
    private double readCountInHomAltSites;
    private double refCountInHomAltSites;
    private double otherAltCountInHomAltSites;
    private double expectedRefInHomAltPerUnitContamination;
    private double expectedRefExcessInHetPerUnitContamination;

    public double getHomAltCount() {
        return homAltCount;
    }
    public double getHetCount() {
        return hetCount;
    }
    public double getExpectedHomAltCount() {
        return expectedHomAltCount;
    }
    public double getExpectedHetCount() {
        return expectedHetCount;
    }


    public void increment(final BiallelicGenotypes.Posterior posteriors, final PileupSummary ps) {
        final double homAltResponsibility = posteriors.get(BiallelicGenotypes.HOM_ALT);
        final double hetResponsibility = posteriors.get(BiallelicGenotypes.HET);

        homAltCount += homAltResponsibility;
        hetCount += hetResponsibility;

        final double homAltPrior = MathUtils.square(ps.getAlleleFrequency());
        final double hetPrior = 2 * ps.getAlleleFrequency() * ( 1 - ps.getAlleleFrequency());
        expectedHomAltCount += homAltPrior;
        expectedHetCount += hetPrior;

        readCountInHomAltSites += homAltResponsibility * ps.getTotalCount();
        refCountInHomAltSites += homAltResponsibility * ps.getRefCount();
        otherAltCountInHomAltSites += homAltResponsibility * ps.getOtherAltCount();

        expectedRefInHomAltPerUnitContamination += homAltResponsibility * ps.getTotalCount() * ps.getRefFrequency();
        expectedRefExcessInHetPerUnitContamination += hetResponsibility * ps.getTotalCount() * ( ps.getRefFrequency() - ps.getAlleleFrequency());
    }

    public void increment(final ContaminationStats other) {
        this.homAltCount += other.homAltCount;
        this.hetCount += other.hetCount;

        this.expectedHomAltCount = other.expectedHomAltCount;
        this.expectedHetCount += other.expectedHetCount;

        this.readCountInHomAltSites += other.readCountInHomAltSites;
        this.refCountInHomAltSites += other.refCountInHomAltSites;
        this.otherAltCountInHomAltSites += other.otherAltCountInHomAltSites;

        this.expectedRefInHomAltPerUnitContamination += other.expectedRefInHomAltPerUnitContamination;
        this.expectedRefExcessInHetPerUnitContamination += other.expectedRefExcessInHetPerUnitContamination;
    }

    public double contaminationEstimate() {
        // if ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final double refInHomAltDueToError = otherAltCountInHomAltSites / 2;
        final double refCountInHomAltDueToContamination = Math.max(refCountInHomAltSites - refInHomAltDueToError, 0);
        return refCountInHomAltDueToContamination / expectedRefInHomAltPerUnitContamination;
    }

    public double standardErrorOfContaminationEstimate() {
        return Math.sqrt(contaminationEstimate() / expectedRefInHomAltPerUnitContamination);
    }

    public static ContaminationStats getStats(final Collection<PileupSummary> pileupSummaries, final double contamination) {
        final ContaminationStats result = new ContaminationStats();
        pileupSummaries.forEach(ps -> result.increment(BiallelicGenotypes.getPosteriors(ps, contamination), ps));
        return result;
    }

    public static ContaminationStats getStats(final List<PileupSummary> pileupSummaries, List<BiallelicGenotypes.Posterior> posteriors) {
        Utils.validateArg(pileupSummaries.size() == posteriors.size(), "Must have one posterior per pileup summary.");
        final ContaminationStats result = new ContaminationStats();
        IntStream.range(0, pileupSummaries.size()).forEach(n -> result.increment(posteriors.get(n), pileupSummaries.get(n)));
        return result;
    }
}
