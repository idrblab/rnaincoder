package edu.put.ma.model.input;

import java.io.File;

import lombok.Getter;

import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableMap;

import edu.put.ma.CommandEnum;
import edu.put.ma.model.services.StructuresSet;
import edu.put.ma.utils.ArrayUtils;
import edu.put.ma.utils.PdbReader;
import edu.put.ma.utils.PreconditionUtils;

public class ComparisonInputModelImpl extends AnalysisInputModelImpl {

    private static final Logger LOGGER = LoggerFactory.getLogger(ComparisonInputModelImpl.class);

    private String referenceStructurePdbFilePath;

    public ComparisonInputModelImpl(final String[] args, final String artifactId) {
        super(args, artifactId);
        secureInitState(artifactId);
    }

    protected ComparisonInputModelImpl(final Builder comparisonInputModelBuilder) {
        super(comparisonInputModelBuilder);
        this.referenceStructurePdbFilePath = comparisonInputModelBuilder.referenceStructurePdbFilePath;
        this.initOptionsMapping();
    }

    @Override
    public boolean isInputInitializedProperly() {
        return super.isInputInitializedProperly() && isCommandLineHasOption("r") && isValid();
    }

    @Override
    public Options constructSpecificOptions() {
        final Options options = super.constructSpecificOptions();
        options.addOption("r", "reference-structure-file-path", true,
                "PDB file path of the reference structure");
        return options;
    }

    @Override
    public boolean isValid() {
        return isAppropriateCommandUsed() && isReferenceStructurePdbFilePathValid();
    }

    @Override
    public StructuresSet constructStructuresSet() {
        final StructuresSet structuresSet = super.constructStructuresSet();
        PdbReader.readSingleStructure(referenceStructurePdbFilePath, structuresSet, true);
        return structuresSet;
    }

    @Override
    public String getInputModelString() {
        final StringBuilder inputModelStringBuilder = new StringBuilder(super.getInputModelString());
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder,
                referenceStructurePdbFilePath, "PDB file path of the reference structure");
        return inputModelStringBuilder.toString();
    }

    @Override
    protected void initOptionsMapping() {
        super.initOptionsMapping();
        optionsMapping.putAll(ImmutableMap.of("referenceStructurePdbFilePath", "-r"));
    }

    @Override
    protected void initState() {
        super.initState();
        setReferenceStructurePdbFilePath();
    }

    @Getter
    public static class Builder extends AnalysisInputModelImpl.Builder {

        private String referenceStructurePdbFilePath;

        public Builder referenceStructurePdbFilePath(final String referenceStructurePdbFilePath) {
            this.referenceStructurePdbFilePath = referenceStructurePdbFilePath;
            return this;
        }

        public AnalysisInputModel build() {
            return new ComparisonInputModelImpl(this);
        }
    }

    private boolean isReferenceStructurePdbFilePathValid() {
        if (StringUtils.isBlank(referenceStructurePdbFilePath)) {
            LOGGER.info("PDB file path of the reference structure should be defined (using '-r' option).");
            return false;
        }
        final File referenceStructurePdbFile = FileUtils.getFile(referenceStructurePdbFilePath);
        PreconditionUtils.checkIfFileExistsAndIsNotADirectory(referenceStructurePdbFile,
                "PDB file of the reference structure");
        return StringUtils.endsWithIgnoreCase(referenceStructurePdbFilePath, ".pdb");
    }

    private boolean isAppropriateCommandUsed() {
        return CommandEnum.Measure.getReferenceStructureBasedMeasures().contains(getCommand().toString())
                || ArrayUtils.getEnumNames(CommandEnum.Sequence.class).contains(getCommand().toString())
                || ArrayUtils.getEnumNames(CommandEnum.Structure3d.class).contains(getCommand().toString());
    }

    private void setReferenceStructurePdbFilePath() {
        referenceStructurePdbFilePath = getOptionString("r");
    }

}
