package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given pileup data from {@link GetPileupSummaries}, calculates the fraction of reads coming from cross-sample contamination.
 *
 * <p>
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 *
 * <p>This tool and GetPileupSummaries together replace GATK3's ContEst.</p>
 *
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculate the fraction of reads coming from cross-sample contamination",
        oneLineSummary = "Calculate the fraction of reads coming from cross-sample contamination",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    public static final String MATCHED_PILEUP_SUMMARIES_TABLE_NAME = "matched";

    private static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    private static final double INITIAL_CONTAMINATION_GUESS = 0.1;

    private static final double CONTAMINATION_CONVERGENCE_THRESHOLD = 0.0001;
    private static final int MAX_ITERATIONS = 10;

    private static final double LOH_LOG_POSTERIOR_THRESHOLD = Math.log(0.1);

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table from the sample in question.")
    private File inputTable;

    @Argument(fullName = MATCHED_PILEUP_SUMMARIES_TABLE_NAME,
            shortName = MATCHED_PILEUP_SUMMARIES_TABLE_NAME,
            doc="An optional input table from a matched sample, presumably with less CNV and less contamination.",
            optional = true)
    private File matchedTable = null;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table")
    private final File outputTable = null;

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = PileupSummary.readPileupSummaries(inputTable);
        final Optional<List<PileupSummary>> matchedSummaries = matchedTable == null ? Optional.empty() :
                Optional.of(PileupSummary.readPileupSummaries(matchedTable));

        if (matchedSummaries.isPresent()) {
            final List<PileupSummary> matched = matchedSummaries.get();
            Utils.validateArg(pileupSummaries.size() == matched.size(), "Matched pileup table is not the same length as sample pileup table");
            final boolean sameSites = IntStream.range(0, pileupSummaries.size())
                    .allMatch(n -> pileupSummaries.get(n).getStart() == matched.get(n).getStart()
                            && pileupSummaries.get(n).getContig().equals(matched.get(n).getContig()));
            Utils.validateArg(sameSites, "Matched pileup table has different sites from sample table.");
        }

        // get non-LoH sites and genotypes from the matched normal if available
        final Pair<int[], List<BiallelicGenotypes.Posterior>> nonLoHIndicesAndPosteriors =
                findNonLoHIndicesAndPosteriors(matchedSummaries.orElseGet(() -> pileupSummaries));

        final List<PileupSummary> nonLoHSites = Arrays.stream(nonLoHIndicesAndPosteriors.getFirst())
                .mapToObj(pileupSummaries::get).collect(Collectors.toList());

        final ContaminationStats nonLohStats = ContaminationStats.getStats(nonLoHSites, nonLoHIndicesAndPosteriors.getSecond());

        ContaminationRecord.writeContaminationTable(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), nonLohStats.contaminationEstimate(), nonLohStats.standardErrorOfContaminationEstimate())), outputTable);

        return "SUCCESS";
    }

    /**
     * Use an HMM to determine which sites are definitely not loss of heterozygosity and what their probabilities of being hom alt are
     * @param pileupSummaries a List of pileup summary sites, some of which might exhibit loss of heterozygosity
     */
    private Pair<int[], List<BiallelicGenotypes.Posterior>> findNonLoHIndicesAndPosteriors(List<PileupSummary> pileupSummaries) {
        logger.info("Looking for non loss of heterozygosity sites in the input pileups.");
        double contamination = INITIAL_CONTAMINATION_GUESS;
        int iteration = 0;
        boolean converged = false;

        while (!converged && iteration++ < MAX_ITERATIONS) {  // loop over contamination convergence
            logger.info("Beginning iteration " + iteration + ".");
            final ContaminationHMM hmm = new ContaminationHMM(contamination);
            final ForwardBackwardAlgorithm.Result<PileupSummary, PileupSummary, Boolean> fbResult =
                    ForwardBackwardAlgorithm.apply(pileupSummaries, pileupSummaries, hmm);
            final double[] lohPosteriors = IntStream.range(0, pileupSummaries.size())
                    .mapToDouble(n -> fbResult.logProbability(n, Boolean.TRUE)).toArray();

            final List<PileupSummary> lohData = IntStream.range(0, pileupSummaries.size())
                    .filter(n -> lohPosteriors[n] > LOH_LOG_POSTERIOR_THRESHOLD)
                    .mapToObj(pileupSummaries::get)
                    .collect(Collectors.toList());

            final int[] nonLohIndices = IntStream.range(0, pileupSummaries.size())
                    .filter(n -> lohPosteriors[n] <= LOH_LOG_POSTERIOR_THRESHOLD).toArray();

            final List<PileupSummary> nonLohData = Arrays.stream(nonLohIndices)
                    .mapToObj(pileupSummaries::get).collect(Collectors.toList());

            final ContaminationStats lohStats = ContaminationStats.getStats(lohData, contamination);
            final ContaminationStats nonLohStats = ContaminationStats.getStats(nonLohData, contamination);

            logger.info("Removing " + lohData.size() + " sites due to possible loss of heterozygosity.");
            logStats(lohStats);

            logger.info("This leaves " + nonLohData.size() + " non-LoH sites.");
            logStats(nonLohStats);

            final double newContamination = nonLohStats.contaminationEstimate();
            converged = Math.abs(contamination - newContamination) < CONTAMINATION_CONVERGENCE_THRESHOLD;
            contamination = newContamination;
            if (converged || iteration == MAX_ITERATIONS) {
                final List<BiallelicGenotypes.Posterior> genotypePosteriors = nonLohData.stream()
                        .map(ps -> BiallelicGenotypes.getPosteriors(ps, newContamination))
                        .collect(Collectors.toList());
                return Pair.create(nonLohIndices, genotypePosteriors);
            }
        }
        throw new GATKException.ShouldNeverReachHereException("This method should have returned within the above loop.");
    }

    private void logStats(ContaminationStats stats) {
        logger.info(String.format("These sites have %.2f hets and %.2f hom alts compared to %.2f and %.2f expected.",
                stats.getHetCount(), stats.getHomAltCount(), stats.getExpectedHetCount(), stats.getExpectedHomAltCount()));
        logger.info(String.format("The contamination estimate from these sites is %.4f.", stats.contaminationEstimate()));
    }
}
