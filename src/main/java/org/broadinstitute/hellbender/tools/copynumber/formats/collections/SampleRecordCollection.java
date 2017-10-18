package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.Sample;
import org.broadinstitute.hellbender.tools.copynumber.formats.SampleMetadataReader;
import org.broadinstitute.hellbender.tools.copynumber.formats.SampleMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents an immutable collection of records associated with a sample, a set of
 * mandatory column headers given by a {@link TableColumnCollection}, and lambdas for
 * reading and writing records.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class SampleRecordCollection<T> implements Sample {
    public static final String COMMENT_PREFIX = "@";

    private final SampleMetadata sampleMetadata;
    private final ImmutableList<T> records;
    private final TableColumnCollection mandatoryColumns;
    private final Function<DataLine, T> dataLineToRecordFunction;
    private final BiConsumer<T, DataLine> recordAndDataLineBiConsumer;

    /**
     * Constructor given the sample metadata, the list of records, the mandatory column headers,
     * and the lambdas for reading and writing records.
     */
    SampleRecordCollection(final SampleMetadata sampleMetadata,
                           final List<T> records,
                           final TableColumnCollection mandatoryColumns,
                           final Function<DataLine, T> dataLineToRecordFunction,
                           final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        this.sampleMetadata = Utils.nonNull(sampleMetadata);
        this.records = ImmutableList.copyOf(Utils.nonNull(records));
        Utils.nonNull(mandatoryColumns);
        Utils.nonEmpty(mandatoryColumns.names());
        this.mandatoryColumns = mandatoryColumns;
        this.dataLineToRecordFunction = Utils.nonNull(dataLineToRecordFunction);
        this.recordAndDataLineBiConsumer = Utils.nonNull(recordAndDataLineBiConsumer);
    }

    /**
     * Constructor given an input file, the mandatory column headers, and the lambdas for reading and writing records.
     * The sample metadata is read from comment strings in the file and the list of records is read using
     * the column headers and the appropriate lambda.
     */
    SampleRecordCollection(final File inputFile,
                           final TableColumnCollection mandatoryColumns,
                           final Function<DataLine, T> dataLineToRecordFunction,
                           final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        IOUtils.canReadFile(inputFile);
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.dataLineToRecordFunction = Utils.nonNull(dataLineToRecordFunction);
        this.recordAndDataLineBiConsumer = Utils.nonNull(recordAndDataLineBiConsumer);

        try (final SampleTSVReader reader = new SampleTSVReader(inputFile)) {
            TableUtils.checkMandatoryColumns(reader.columns(), mandatoryColumns, UserException.BadInput::new);
            sampleMetadata = reader.readSampleMetadata();
            records = reader.stream().collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf));
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    @Override
    public final SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    public final int size() {
        return records.size();
    }

    /**
     * @return  an immutable view of the records contained in the collection
     */
    public final List<T> getRecords() {
        return records;
    }

    /**
     * Writes the sample name in a comment string and the records to file.
     */
    public void write(final File outputFile) {
        try (final SampleTSVWriter writer = new SampleTSVWriter(outputFile, sampleMetadata)) {
            writer.writeSampleMetadata();
            writer.writeAllRecords(records);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final SampleRecordCollection<?> that = (SampleRecordCollection<?>) o;
        return sampleMetadata.equals(that.sampleMetadata) && records.equals(that.records);
    }

    @Override
    public int hashCode() {
        int result = sampleMetadata.hashCode();
        result = 31 * result + records.hashCode();
        return result;
    }

    private final class SampleTSVReader extends TableReader<T> {
        private final File file;

        private SampleTSVReader(final File file) throws IOException {
            super(file);
            this.file = file;
        }

        private SampleMetadata readSampleMetadata() {
            return new SampleMetadata(new SAMFileHeader());
        }

        @Override
        protected T createRecord(final DataLine dataLine) {
            Utils.nonNull(dataLine);
            Utils.validateArg(dataLine.get(0).startsWith(COMMENT_PREFIX),
                    "The value of the first column cannot start with the comment prefix " +COMMENT_PREFIX);
            return dataLineToRecordFunction.apply(dataLine);
        }
    }

    private final class SampleTSVWriter extends TableWriter<T> {
        private final SampleMetadata sampleMetadata;

        SampleTSVWriter(final File file,
                        final SampleMetadata sampleMetadata) throws IOException {
            super(file, mandatoryColumns, COMMENT_PREFIX);
            this.sampleMetadata = Utils.nonNull(sampleMetadata);
        }

        void writeSampleMetadata() {
            try {
                writeComment();
            } catch (final IOException e) {
                throw new UserException("Could not write sample name.");
            }
        }

        @Override
        protected void composeLine(final T record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);
            recordAndDataLineBiConsumer.accept(record, dataLine);
        }
    }
}