package org.broadinstitute.hellbender.tools.copynumber.temporary;

import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;

public class HDF5ReadCountCollection {
    private static final String SAMPLE_NAME_PATH = "/sample_name/value";
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String READ_COUNTS_PATH = "/read_counts/values";

    private final HDF5File file;
    private Lazy<List<SimpleInterval>> intervals;

    public HDF5ReadCountCollection(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
    }

    public String getSampleName() {
        return file.readStringArray(SAMPLE_NAME_PATH)[0];
    }

    public List<SimpleInterval> getIntervals() {
        return intervals.get();
    }

    public RealMatrix getReadCounts() {
        return new Array2DRowRealMatrix(file.readDoubleMatrix(READ_COUNTS_PATH));
    }

    public static void write(final File outFile,
                             final String sampleName,
                             final List<SimpleInterval> intervals,
                             final double[][] values) {
        Utils.nonNull(outFile);
        Utils.nonNull(sampleName);
        Utils.nonNull(intervals);
        Utils.nonNull(values);

        if (values[0].length != intervals.size()) {
            throw new GATKException("The shape of the values array (" + values.length + " x " + values[0].length +
                    ") does not match the number of intervals (" + intervals.size() + ").");
        }

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5ReadCountCollection hdf5ReadCountCollection = new HDF5ReadCountCollection(file);
            hdf5ReadCountCollection.writeName(SAMPLE_NAME_PATH, sampleName);
            hdf5ReadCountCollection.writeIntervals(intervals);
            hdf5ReadCountCollection.writeReadCounts(values);
        }
    }

    private <T extends SimpleInterval> void writeIntervals(final List<T> intervals) {
        HDF5Utils.writeIntervals(file, INTERVALS_GROUP_NAME, intervals);
    }

    private void writeName(final String path, final String name) {
        file.makeStringArray(path, name);
    }

    private void writeReadCounts(final double[][] readCounts) {
        file.makeDoubleMatrix(READ_COUNTS_PATH, readCounts);
    }
}
