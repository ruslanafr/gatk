package org.broadinstitute.hellbender.utils.realignmentfilter;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class Realigner {
    final private BwaMemAligner aligner;

    public Realigner(final RealignmentFilterArgumentCollection rfac) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        aligner = new BwaMemAligner(index);


    }
}
