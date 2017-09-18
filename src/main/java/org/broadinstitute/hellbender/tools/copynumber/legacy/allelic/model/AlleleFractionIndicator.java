package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

/**
 * Possible hidden states of sites in allele fraction model.
 *
 * NOTE: these are not part of AlleleFractionState because AlleleFractionState pertains to a collapsed model in which latent variables
 * (this indicator and allelic biases) have been marginalized out.  These hidden states are, however, necessary
 * for initializing the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public enum AlleleFractionIndicator {
    HOM_ALT, HOM_REF, HET_ALT_MINOR, HET_REF_MINOR, OUTLIER
}
