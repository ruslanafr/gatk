package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.*;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.VariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.samples.Sample;

import java.io.File;
import java.util.*;

/**
 * Annotate variant calls with context information
 *
 * <p>
 * This tool is designed to annotate variant calls based on their context (as opposed to functional annotation).
 * Various annotation modules are available; see the "Annotation Modules" page linked in the Tool Documentation sidebar for a complete list.
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to annotate and optionally one or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * An annotated VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <br />
 *
 * <h4>Annotate a VCF with dbSNP IDs and depth of coverage for each sample</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T VariantAnnotator \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -A Coverage \
 *   -L input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>Annotate a VCF with allele frequency by an external resource. Annotation will only occur if there is allele concordance between the resource and the input VCF </h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T VariantAnnotator \
 *   -I input.bam \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -L input.vcf \
 *   --resource:foo resource.vcf \
 *   -E foo.AF \
 *   --resourceAlleleConcordance
 * </pre>
 *
 * <h4>Annotate with AF and FILTER fields from an external resource </h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T VariantAnnotator \
 *   -V input.vcf \
 *   -o output.vcf \
 *   --resource:foo resource.vcf \
 *   --expression foo.AF \
 *   --expression foo.FILTER
 * </pre>
 *
 */
public class VariantAnnotator extends VariantWalker {
    private VariantContextWriter vcfWriter;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName = "comp", shortName = "comp", doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = Collections.emptyList();

    /**
     * An external resource VCF file or files from which to annotate.
     *
     * Use this option to add annotations from a resource file to the output.
     * For example, if you want to annotate your callset with the AC field value from a VCF file named
     * 'resource_file.vcf', you tag it with '-resource:my_resource resource_file.vcf' and you additionally specify
     * '-E my_resource.AC' (-E is short for --expression, also documented on this page). In the resulting output
     * VCF, any records for which there is a record at the same position in the resource file will be annotated with
     * 'my_resource.AC=N'. Note that if there are multiple records in the resource file that overlap the given
     * position, one is chosen randomly. Check for allele concordance if using --resourceAlleleConcordance, otherwise
     * the match is based on position only.
     */
    @Input(fullName="resource", shortName = "resource", doc="External resource VCF file", required=false)
    public List<RodBinding<VariantContext>> resources = Collections.emptyList();
    public List<RodBinding<VariantContext>> getResourceRodBindings() { return resources; }

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The file to whcih variants should be written", optional=false)
    private File outputFile;

    /**
     * See the --list argument to view available annotations.
     */
    @ArgumentCollection
    public VariantAnnotationArgumentCollection variantAnnotationArgumentCollection = new VariantAnnotationArgumentCollection(
            Arrays.asList(StandardAnnotation.class.getSimpleName()),
            Collections.emptyList(),
            Collections.emptyList());

    /**
     * This option enables you to add annotations from one VCF to another.
     *
     * For example, if you want to annotate your callset with the AC field value from a VCF file named
     * 'resource_file.vcf', you tag it with '-resource:my_resource resource_file.vcf' (see the -resource argument, also
     * documented on this page) and you specify '-E my_resource.AC'. In the resulting output VCF, any records for
     * which there is a record at the same position in the resource file will be annotated with 'my_resource.AC=N'.
     * INFO field data, ID, ALT, and FILTER fields may be used as expression values.
     * Note that if there are multiple records in the resource file that overlap the given position, one is chosen
     * randomly.
     */
    @Argument(fullName="expression", shortName="E", doc="One or more specific expressions to apply to variant calls", optional=true)
    protected Set<String> expressionsToUse = new ObjectOpenHashSet();

    /**
     * If this argument is specified, add annotations (specified by --expression) from an external resource
     * (specified by --resource) to the input VCF (specified by --variant) only if the alleles are
     * concordant between input and the resource VCFs. Otherwise, always add the annotations.
     */
    @Argument(fullName="resourceAlleleConcordance", shortName="rac", doc="Check for allele concordances when using an external resource VCF file", optional=true)
    protected Boolean expressionAlleleConcordance = false;

    /**
     * You can use the -XL argument in combination with this one to exclude specific annotations.Note that some
     * annotations may not be actually applied if they are not applicable to the data provided or if they are
     * unavailable to the tool (e.g. there are several annotations that are currently not hooked up to
     * HaplotypeCaller). At present no error or warning message will be provided, the annotation will simply be
     * skipped silently. You can check the output VCF header to see which annotations were actually applied (although
     * this does not guarantee that the annotation was applied to all records in the VCF, since some annotations have
     * additional requirements, e.g. minimum number of samples or heterozygous sites only -- see the documentation
     * for individual annotations' requirements).
     */
    @Argument(fullName="useAllAnnotations", shortName="all", doc="Use all possible annotations (not for the faint of heart)", optional=true)
    protected Boolean USE_ALL_ANNOTATIONS = false;

    /**
     * Note that the --list argument requires a fully resolved and correct command-line to work. As an alternative,
     * you can use ListAnnotations (see Help Utilities).
     */
    @Argument(fullName="list", shortName="ls", doc="List the available annotations and exit", optional=true)
    protected Boolean LIST = false;

    /**
     * By default, a dbSNP ID is added only when the ID field in the variant record is empty (not already annotated).
     * This argument allows you to override that behavior, and appends the new ID to the existing one. This is used
     * in conjunction with the -dbsnp argument.
     */
    @Argument(fullName="alwaysAppendDbsnpId", shortName="alwaysAppendDbsnpId", doc="Add dbSNP ID even if one is already present", optional=true)
    protected Boolean ALWAYS_APPEND_DBSNP_ID = false;
    public boolean alwaysAppendDbsnpId() { return ALWAYS_APPEND_DBSNP_ID; }

    /**
     * The genotype quality (GQ) threshold above which the mendelian violation ratio should be annotated.
     */
    @Argument(fullName="MendelViolationGenotypeQualityThreshold",shortName="mvq",optional=true,doc="GQ threshold for annotating MV ratio")
    public double minGenotypeQualityP = 0.0;

    private VariantAnnotatorEngine annotatorEngine;

    /**
     * Prepare the output file and the list of available features.
     */
    public void onTraversalStart() {
        //TODO figure out how to do this with barclay?
//        if ( LIST ) {
//            AnnotationHelpUtils.listAnnotations();
//            System.exit(0);
//        }

        // get the list of all sample names from the variant VCF input rod, if applicable
        final  List<String> samples = getHeaderForVariants().getSampleNamesInOrder();

        if ( USE_ALL_ANNOTATIONS ) {
            annotatorEngine = VariantAnnotatorEngine.ofAllMinusExcluded(variantAnnotationArgumentCollection.annotationsToExclude, dbsnp.dbsnp, comps);
        } else {
            annotatorEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(variantAnnotationArgumentCollection, dbsnp.dbsnp, comps);
        }
        //TODO add expressions?
//        annotatorEngine.initializeExpressions(expressionsToUse);...
//        annotatorEngine.setExpressionAlleleConcordance(expressionAlleleConcordance);

        // setup the header fields
        // note that if any of the definitions conflict with our new ones, then we want to overwrite the old ones
        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.addAll(annotatorEngine.getVCFAnnotationDescriptions(false)); //TODO allele specific annotations
        for ( final VCFHeaderLine line : getHeaderForVariants().getLin ) {
            if ( isUniqueHeaderLine(line, hInfo) )
                hInfo.add(line);
        }
        // for the expressions, pull the info header line from the header of the resource rod
        for ( final VariantAnnotatorEngine.VAExpression expression : annotatorEngine.getRequestedExpressions() ) {
            // special case the ID field
            if ( expression.fieldName.equals("ID") ) {
                hInfo.add(new VCFInfoHeaderLine(expression.fullName, 1, VCFHeaderLineType.String, "ID field transferred from external VCF resource"));
                continue;
            }
            VCFInfoHeaderLine targetHeaderLine = null;
            for ( final VCFHeaderLine line : GATKVCFUtils.getHeaderFields(getToolkit(), Arrays.asList(expression.binding.getName())) ) {
                if ( line instanceof VCFInfoHeaderLine ) {
                    final VCFInfoHeaderLine infoline = (VCFInfoHeaderLine)line;
                    if ( infoline.getID().equals(expression.fieldName) ) {
                        targetHeaderLine = infoline;
                        break;
                    }
                }
            }

            if ( targetHeaderLine != null ) {
                if ( targetHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER )
                    hInfo.add(new VCFInfoHeaderLine(expression.fullName, targetHeaderLine.getCount(), targetHeaderLine.getType(), targetHeaderLine.getDescription()));
                else
                    hInfo.add(new VCFInfoHeaderLine(expression.fullName, targetHeaderLine.getCountType(), targetHeaderLine.getType(), targetHeaderLine.getDescription()));
            } else {
                hInfo.add(new VCFInfoHeaderLine(expression.fullName, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Value transferred from another external VCF resource"));
            }
        }

        //annotatorEngine.makeHeaderInfoMap(hInfo);
       // annotatorEngine.invokeAnnotationInitializationMethods(hInfo);

        VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(vcfHeader);
    }

    public static boolean isUniqueHeaderLine(VCFHeaderLine line, Set<VCFHeaderLine> currentSet) {
        if ( !(line instanceof VCFCompoundHeaderLine) )
            return true;

        for ( VCFHeaderLine hLine : currentSet ) {
            if ( hLine instanceof VCFCompoundHeaderLine && ((VCFCompoundHeaderLine)line).sameLineTypeAndName((VCFCompoundHeaderLine)hLine) )
                return false;
        }

        return true;
    }

    /**
     * We want reads that span deletions
     *
     * @return true
     */
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    /**
     * For each site of interest, annotate based on the requested annotation types
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {

        // if the reference base is not ambiguous, we can annotate
        Map<String, AlignmentContext> stratifiedContexts;
        if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(readsContext.());
        }
        if ( BaseUtils.simpleBaseToBaseIndex(refContext.getBase()) != -1 ) {
            List<Sample> samplesList = vc.getSampleNamesOrderedByName();
            ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>( vc.getSampleNamesOrderedByName(), vc.getAlleles(),
                    AssemblyBasedCallerUtils.splitReadsBySample(samplesList, readsHeader, readsContext);)
            VariantContext annotatedVC = annotatorEngine.annotateContext(vc, fc, refContext, stratifiedContexts, a -> true);
            vcfWriter.add(annotatedVC);
        }
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     */
    @Override
    public Object onTraversalSuccess() {
//        logger.info("Processed " + result + " loci.\n");
        return null;
    }
}

