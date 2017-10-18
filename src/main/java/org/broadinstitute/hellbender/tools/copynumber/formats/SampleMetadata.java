package org.broadinstitute.hellbender.tools.copynumber.formats;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class SampleMetadata {
    private final String sampleName;
    private final SAMSequenceDictionary sequenceDictionary;

    public SampleMetadata(final SAMFileHeader header) {
        Utils.nonNull(header);
        Utils.nonNull(header.getReadGroups());
        final List<String> sampleNames = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());
        Utils.validateArg(sampleNames.size() == 1,"Must create metadata from data with exactly one sample name.");
        this.sampleName = sampleNames.get(0);
        Utils.validateArg(!header.getSequenceDictionary().isEmpty(), "Must create metadata from data with a non-empty sequence dictionary.");
        this.sequenceDictionary = header.getSequenceDictionary();
    }

    public String getSampleName() {
        return sampleName;
    }

    public SAMSequenceDictionary getSequenceDictionary() {
        return new SAMSequenceDictionary(sequenceDictionary.getSequences());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SampleMetadata that = (SampleMetadata) o;
        return sampleName.equals(that.sampleName) && sequenceDictionary.isSameDictionary(that.sequenceDictionary);
    }

    @Override
    public int hashCode() {
        int result = sampleName.hashCode();
        result = 31 * result + sequenceDictionary.hashCode();
        return result;
    }
}
