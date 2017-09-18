package org.broadinstitute.hellbender.tools.copynumber.temporary;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class HDF5ReadCountCollectionUnitTest extends BaseTest {

    @Test
    public void basicTest() {
        final File outputFile = createTempFile("HDF5ReadCountCollection", ".hdf5");
        final String sampleName = "SAMPLE1";
        final List<SimpleInterval> newIntervals = new ArrayList<>();
        newIntervals.add(new SimpleInterval("1", 1000, 2000));
        newIntervals.add(new SimpleInterval("1", 5000, 6000));
        final double[][] newValues = {{2, 10}};

        // The output file already exists at this point, since it is a temp file.
        HDF5ReadCountCollection.write(outputFile, sampleName, newIntervals, newValues);

        final HDF5ReadCountCollection rcc = new HDF5ReadCountCollection(new HDF5File(outputFile));

        Assert.assertEquals(rcc.getSampleName(), sampleName);

        Assert.assertEquals(rcc.getIntervals(), newIntervals);
        Assert.assertFalse(rcc.getIntervals() == newIntervals);

        Assert.assertTrue(IntStream.range(0, newValues.length).allMatch(i -> newValues[i][0] == rcc.getReadCounts().getData()[0][i]));

        Assert.assertEquals(rcc.getReadCounts().transpose().getData().length, newValues.length);
        Assert.assertEquals(rcc.getReadCounts().transpose().getData()[0].length, newValues[0].length);
        Assert.assertFalse(rcc.getReadCounts().transpose().getData() == newValues);
    }
}
