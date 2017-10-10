package org.broadinstitute.hellbender.utils.realignmentfilter;

import org.broadinstitute.barclay.argparser.Argument;

public class RealignmentFilterArgumentCollection {

    /**
     * BWA-mem index image.
     */
    @Argument(fullName = "bwa-mem-index-image", shortName = "index", doc = "BWA-mem index image")
    public String bwaMemIndexImage;

}
