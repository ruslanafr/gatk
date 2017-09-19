package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegment;

import java.util.List;

public class CalledCopyRatioSegment extends CopyRatioSegment {
    enum Call {
        AMPLIFICATION_CALL("+"),
        DELETION_CALL("-"),
        NEUTRAL_CALL("0");

        private final String outputString;

        Call(final String outputString) {  this.outputString = outputString; }

        public String getOutputString() {
            return outputString;
        }
    }

    private final Call call;

    public CalledCopyRatioSegment(final CopyRatioSegment segment,
                                  final Call call) {
        super(segment.getInterval(), segment.getNumPoints(), segment.getMeanLog2CopyRatio());
        this.call = call;
    }

    public Call getCall() {
        return call;
    }
}
