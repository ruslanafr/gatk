package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CopyRatioCollection {
    private final String sampleName;
    private final List<CopyRatio> copyRatios;

    /**
     * Constructor from file.
     * @param inputFile file to read from
     * @throws UserException.CouldNotReadInputFile if input file cannot be read
     */
    public CopyRatioCollection(final File inputFile) {
        IOUtils.canReadFile(inputFile);

        try (final CopyRatioReader reader = new CopyRatioReader(inputFile)) {
            sampleName = reader.getSampleName();
            copyRatios = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public CopyRatioCollection(final String sampleName,
                               final List<SimpleInterval> intervals,
                               final RealMatrix copyRatioValues) {
        Utils.nonNull(sampleName);
        Utils.nonNull(intervals);
        Utils.nonNull(copyRatioValues);
        Utils.validateArg(copyRatioValues.getRowDimension() == 1,
                "Copy-ratio values must contain only a single row.");
        Utils.validateArg(intervals.size() == copyRatioValues.getColumnDimension(),
                "Number of intervals and columns in copy-ratio values must match.");
        this.sampleName = sampleName;
        this.copyRatios = IntStream.range(0, intervals.size()).boxed()
                .map(i -> new CopyRatio(intervals.get(i), copyRatioValues.getRow(0)[i]))
                .collect(Collectors.toList());
    }

    public String getSampleName() {
        return sampleName;
    }

    /**
     * Returns an unmodifiable view of the list of {@link CopyRatio}s.
     */
    public List<CopyRatio> getCopyRatios() {
        return Collections.unmodifiableList(copyRatios);
    }

    /**
     * Returns a new array of the copy-ratio values.
     */
    public double[] getCopyRatioValues() {
        return copyRatios.stream().mapToDouble(CopyRatio::getLog2CopyRatioValue).toArray();
    }

    /**
     * Writes copy ratios to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     * @throws UserException.CouldNotCreateOutputFile if output file cannot be created
     */
    public void write(final File outputFile) {
        try (final CopyRatioWriter writer = new CopyRatioWriter(outputFile, sampleName)) {
            writer.writeSampleName();
            writer.writeAllRecords(copyRatios);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }
}
